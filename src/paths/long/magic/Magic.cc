///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "PrintAlignment.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "paths/long/ultra/FounderAlignment.h"
#include "paths/long/ultra/ThreadedBlocks.h"

// QueryLookupTable K=12 MM=12 MC=0.005 SEQS=Jan30.2014.3/all.fastb L=scs.lookup 
//      SMITH_WAT=True FW_ONLY=True PARSEABLE=True > Jan30.2014.3/all.aligns

template<int K> void MakeKmerLookup0SingleX( const vecbasevector& unibases,
     vec< triple<kmer<K>,int,int> >& kmers_plus )
{    vec<int64_t> starts;
     starts.push_back(0);
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    const basevector& u = unibases[i];
          starts.push_back( starts.back( ) + Max( 0, u.isize( ) - K + 1 ) );    }
     kmers_plus.resize( starts.back( ) );
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    const basevector& u = unibases[i];
          for ( int j = 0; j <= u.isize( ) - K; j++ )
          {    int64_t r = starts[i] + j;
               kmers_plus[r].first.SetToSubOf( u, j );
               kmers_plus[r].second = i;
               kmers_plus[r].third = j;    }    }
     Sort(kmers_plus);    }

template<int K> void MakeBlocksX( 
     const vecbasevector& reads,
     const vec< vec< pair<int,int> > >& a,    // kmer aligns to read 0
     threaded_blocks& tb, ostream& out )                    // output
{
     // Heuristics.

     const int max_shift = 5;
     // ****************************************************************************
     const int min_orbit = 2; // WAS 2 *********************************************
     const double win_ratio = 4.0;
     // ****************************************************************************
     const int min_block = 5; // WAS 50 ********************************************

     // Test input requirements.

     for ( int id = 1; id < (int) reads.size( ); id++ )
     {    ForceAssert( a[id].nonempty( ) );
          for ( int j = 0; j < a[id].isize( ) - 1; j++ )
          {    ForceAssertLt( a[id][j].first, a[id][j+1].first );
               ForceAssertLt( a[id][j].second, a[id][j+1].second );    }    }

     // Define pointed kmers.

     int N = reads.size( );
     vec< pair<kmer<K>,int> > pkmers0;
     for ( int p = 0; p <= reads[0].isize( ) - K; p++ )
     {    kmer<K> x;
          x.SetToSubOf( reads[0], p );
          pkmers0.push( x, p );    }
     for ( int id = 1; id < N; id++ )
     {    for ( int p2 = 0; p2 < a[id][0].first; p2++ )
          {    int p1 = a[id][0].second - a[id][0].first + p2;
               if ( p1 < 0 || p1 >= reads[0].isize( ) ) continue;
               kmer<K> x;
               x.SetToSubOf( reads[id], p2 );
               pkmers0.push( x, p1 );    }
          for ( int j = 0; j < a[id].isize( ); j++ )
          {    int start2 = a[id][j].first, stop2;
               if ( j < a[id].isize( ) - 1 ) stop2 = a[id][j+1].first;
               else stop2 = reads[id].size( );
               for ( int p2 = start2; p2 < stop2; p2++ )
               {    int p1 = a[id][j].second + p2 - start2;
                    if ( p1 < 0 || p1 >= reads[0].isize( ) ) continue;
                    if ( p2 > reads[id].isize( ) - K ) continue;
                    kmer<K> x;
                    x.SetToSubOf( reads[id], p2 );
                    pkmers0.push( x, p1 );    }    }    }
     Sort(pkmers0);

     // Define an equivalence relation on the pointed kmers, by joining (x,p) to
     // (x,p') if |p-p'| <= 5.

     equiv_rel e( pkmers0.size( ) );
     for ( int j = 1; j < pkmers0.isize( ); j++ )
     {    if ( pkmers0[j].first == pkmers0[j-1].first
               && pkmers0[j].second - pkmers0[j-1].second <= max_shift )
          {    e.Join( j-1, j );    }    }

     // Delete orbits of size < 3.

     vec<int> reps;
     e.OrbitReps(reps);
     vec<Bool> small_orbit( reps.size( ), False );
     for ( int j = 0; j < reps.isize( ); j++ )
     {    vec<int> o;
          e.Orbit( reps[j], o );
          if ( o.isize( ) < min_orbit ) small_orbit[j] = True;    }
     EraseIf( reps, small_orbit );
     int nreps = reps.size( );
     PRINT_TO( out, nreps ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     // Find connections, make a graph.

     vec< pair<kmer<K>,int> > qkmers0;
     vec<int> rep_id;
     for ( int j = 0; j < reps.isize( ); j++ )
     {    vec<int> o;
          e.Orbit( reps[j], o );
          for ( int l = 0; l < o.isize( ); l++ )
          {    qkmers0.push_back( pkmers0[ o[l] ] );
               rep_id.push_back(j);    }    }
     SortSync( qkmers0, rep_id );
     vec< vec<int> > from(nreps), to(nreps);
     for ( int id = 1; id < N; id++ )
     {    vec< pair<kmer<K>,int> > tkmers0;
          tkmers0.reserve( reads[id].isize( ) - K + 1 );
          kmer<K> x;
          for ( int j = 0; j < a[id].isize( ); j++ )
          {    int start2 = a[id][j].first, stop2;
               if ( j < a[id].isize( ) - 1 ) stop2 = a[id][j+1].first;
               else stop2 = reads[id].size( );
               for ( int p2 = start2; p2 < stop2; p2++ )
               {    int p1 = a[id][j].second + p2 - start2;
                    if ( p2 > reads[id].isize( ) - K ) continue;
                    x.SetToSubOf( reads[id], p2 );
                    tkmers0.push( x, p1 );    }    }
          for ( int j = 0; j < tkmers0.isize( ) - 1; j++ )
          {    int x1 = BinPosition( qkmers0, tkmers0[j] );
               if ( x1 < 0 ) continue;
               int x2 = BinPosition( qkmers0, tkmers0[j+1] );
               if ( x2 < 0 ) continue;
               int v1 = rep_id[x1], v2 = rep_id[x2];
               from[v1].push_back(v2), to[v2].push_back(v1);    }    }
     for ( int v = 0; v < nreps; v++ )
     {    UniqueSort(from[v]), UniqueSort(to[v]);    }
     digraph G( from, to );

     // Project onto read 0 and remove dominated kmers.

     vec<Bool> bad( G.N( ), False );
     vec<int> start, stop, osize;
     for ( int v = 0; v < nreps; v++ )
     {    vec<int> o, pos;
          e.Orbit( reps[v], o );
          osize.push_back( o.size( ) );
          for ( int l = 0; l < o.isize( ); l++ )
               pos.push_back( pkmers0[ o[l] ].second );
          Sort(pos);
          start.push_back( pos.front( ) );
          stop.push_back( pos.back( ) );    }
     vec<int> ids( nreps, vec<int>::IDENTITY ), max_stop(nreps);
     SortSync( start, stop, ids );
     int ms = 0;
     for ( int j = 0; j < nreps; j++ )
     {    max_stop[j] = ms;
          ms = Max( ms, stop[j] );    }
     for ( int j1 = 0; j1 < nreps; j1++ )
     {    int v1 = ids[j1];
          for ( int j2 = j1 + 1; j2 < nreps; j2++ )
          {    if ( start[j2] > stop[j1] + max_shift ) break;
               int v2 = ids[j2];
               if ( osize[v2] > win_ratio * osize[v1] ) bad[v1] = True;    }
          for ( int j2 = j1 - 1; j2 >= 0; j2-- )
          {    if ( stop[j2] + max_shift < start[j1] ) continue;
               int v2 = ids[j2];
               if ( osize[v2] > win_ratio * osize[v1] ) bad[v1] = True;
               if ( max_stop[j2] + max_shift < start[j1] ) break;    }    }
     PRINT2_TO( out, G.N( ), Sum(bad) ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     // Define maximal unbranched stretches in the graph.

     vec< vec<int> > lines;
     vec<Bool> used( nreps, False );
     for ( int v = 0; v < nreps; v++ )
     {    if ( bad[v] ) continue;
          if ( used[v] ) continue;
          vec<int> line;
          line.push_back(v);
          used[v] = True;
          int w = v;
          while(1)
          {    vec<int> tos;
               for ( int j = 0; j < G.To(w).isize( ); j++ )
               {    int x = G.To(w)[j];
                    if ( !bad[x] ) tos.push_back(x);    }
               if ( !tos.solo( ) ) break;
               w = tos[0];
               if ( Member( line, w ) ) break;
               if ( used[w] ) break; // ???
               int froms = 0;
               for ( int j = 0; j < G.From(w).isize( ); j++ )
               {    int x = G.From(w)[j];
                    if ( !bad[x] ) froms++;    }
               if ( froms != 1 ) break;
               line.push_front(w);
               used[w] = True;    }
          while(1)
          {    vec<int> froms;
               for ( int j = 0; j < G.From(w).isize( ); j++ )
               {    int x = G.From(w)[j];
                    if ( !bad[x] ) froms.push_back(x);    }
               if ( !froms.solo( ) ) break;
               w = froms[0];
               int tos = 0;
               for ( int j = 0; j < G.To(w).isize( ); j++ )
               {    int x = G.To(w)[j];
                    if ( !bad[x] ) tos++;    }
               if ( tos != 1 ) break;
               if ( Member( line, w ) ) break;
               if ( used[w] ) break; // ???
               line.push_back(w);
               used[w] = True;    }
          lines.push_back(line);    }
     PRINT_TO( out, lines.size( ) ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     // Compute and order blocks.

     vec< vec<int> > blocksx;
     vec<int> line_id;
     for ( int m = 0; m < lines.isize( ); m++ )
     {    const vec<int>& L = lines[m];
          if ( L.isize( ) < min_block ) continue;
          blocksx.push_back(L);
          line_id.push_back(m);    }
     vec<int> bpos;
     for ( int m = 0; m < blocksx.isize( ); m++ )
     {    const vec<int>& L = blocksx[m];
          vec<int> pos;
          for ( int k = 0; k < L.isize( ); k++ )
          {    int v = L[k];
               vec<int> o;
               e.Orbit( reps[v], o );
               for ( int l = 0; l < o.isize( ); l++ )
                    pos.push_back( pkmers0[ o[l] ].second );    }
          Sort(pos);
          bpos.push_back( pos[ pos.size( )/2 ] );    }
     SortSync( bpos, blocksx, line_id );

     // Translate pkmers0 into pkmers.

     vec< pair<basevector,int> > pkmers( pkmers0.size( ) );
     for ( int j = 0; j < pkmers0.isize( ); j++ )
     {    pkmers0[j].first.GetBasevector( pkmers[j].first );
          pkmers[j].second = pkmers0[j].second;    }

     // Two successive blocks could overlap by up to K-1.  Trim their ends so that 
     // they don't overlap and in fact so that there is a small gap between them.
     // The small gap is desirable between we don't want threads to have negative
     // lengths.  Of course when could further trim the ends as needed at the point
     // when we create the threads.

     // ****************************************************************************
     /*
     const int target_gap = 3;
     for ( int j = 0; j < blocksx.isize( ) - 1; j++ )
     {    vec<int> &L1 = blocksx[j], &L2 = blocksx[j+1];
          int over;
          const basevector& b1 = pkmers[ reps[ L1.back( ) ] ].first;
          const basevector& b2 = pkmers[ reps[ L2.front( ) ] ].first;
          for ( over = K-1; over >= 1; over-- )
               if ( b1.Overlap( b2, over ) ) break;
          int trim = over + target_gap;
          for ( int j = 0; j < trim; j++ )
          {    if ( L1.size( ) >= L2.size( ) ) L1.resize( L1.isize( ) - 1 );
               else L2.erase( L2.begin( ) );    }    }
     */

     // Compute block positions again.

     vec<ho_interval> bposh;
     for ( int m = 0; m < blocksx.isize( ); m++ )
     {    const vec<int>& L = blocksx[m];
          vec<int> pos;
          for ( int k = 0; k < L.isize( ); k++ )
          {    int v = L[k];
               vec<int> o;
               e.Orbit( reps[v], o );
               for ( int l = 0; l < o.isize( ); l++ )
                    pos.push_back( pkmers[ o[l] ].second );    }
          Sort(pos);
          bposh.push( pos.front( ), pos.back( ) + K );    }

     // Compute the sequences associated to the blocks.

     vec<basevector> blocks;
     for ( int m = 0; m < blocksx.isize( ); m++ )
     {    const vec<int>& L = blocksx[m];
          basevector seq;
          for ( int k = 0; k < L.isize( ); k++ )
          {    int v = L[k];
               vec<int> o;
               e.Orbit( reps[v], o );
               const basevector& b = pkmers[ o[0] ].first;
               if ( seq.size( ) == 0 ) seq = b;
               else 
               {    seq.resize( seq.isize( ) - (K-1) );
                    seq = Cat( seq, b );    }    }
          blocks.push_back(seq);    }
     PRINT_TO( out, blocks.size( ) ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     // Predict initial positions of blocks on each read.

     vec< vec<Bool> > block_defined(N);
     vec< vec<int> > block_start(N), block_stop(N);
     vec< vec<double> > block_err_rate(N);
     for ( int id = 0; id < N; id++ )
     {
          // Create map of positions on read 0 to positions on read id.

          vec<int> to_id( reads[0].size( ) + 1, -1 );
          for ( int p = 0; p < a[id][0].second; p++ )
               to_id[p] = p + a[id][0].first - a[id][0].second;
          for ( int j = 0; j < a[id].isize( ) - 1; j++ )
          {    for ( int p = a[id][j].second; p < a[id][j+1].second; p++ )
                    to_id[p] = p + a[id][j].first - a[id][j].second;    }
          for ( int p = a[id].back( ).second; p <= reads[0].isize( ); p++ )
               to_id[p] = p + a[id].back( ).first - a[id].back( ).second;

          // Now map the blocks.

          block_defined[id].resize( blocks.size( ), False );
          block_start[id].resize( blocks.size( ) );
          block_stop[id].resize( blocks.size( ) );
          block_err_rate[id].resize( blocks.size( ) );
          for ( int m = 0; m < blocksx.isize( ); m++ )
          {    const vec<int>& L = blocksx[m];
               int left = 0, right = 0;
               if ( bposh[m].Start( ) < 0 ) left = -bposh[m].Start( );
               if ( bposh[m].Stop( ) > reads[0].isize( ) )
                    right = bposh[m].Stop( ) - reads[0].isize( );
               int start = left + to_id[ Max( 0, bposh[m].Start( ) ) ];
               int stop = right 
                    + to_id[ Min( reads[0].isize( ), bposh[m].Stop( ) ) ] + K - 1;

               // Align block.  For now assume fully on read id.

               if ( start < 0 || stop > reads[id].isize( ) ) continue;
               if ( !( start < stop ) ) continue;
               align a;
               const int add = 5;
               start = start - add;
               stop = stop + add;
               if ( start < 0 ) start = 0;
               if ( stop > reads[id].isize( ) ) stop = reads[id].size( );
               basevector r( reads[id], start, stop - start );

               int errors;
               const int bw_add = 10 + blocks[m].isize( )/50;
               int offset = -( r.isize( ) - blocks[m].isize( ) ) / 2;
               int bandwidth = Max( 0, ( r.isize( ) - blocks[m].isize( ) ) / 2 );
               int errs = SmithWatBandedA( blocks[m], r, offset, bandwidth + bw_add,
                    a, errors, 0, 1, 1 );

               // SmithWatFreeSym( blocks[m], r, a, false, false, 1, 1 );

               int start2 = start + a.pos2( ), stop2 = start + a.Pos2( );
               block_err_rate[id][m] = double(errs) / double( stop2 - start2 );
               block_defined[id][m] = True;
               block_start[id][m] = start2, block_stop[id][m] = stop2;    }    }
     vec< vec<int> > block_start_orig(block_start), block_stop_orig(block_stop);

     // Delete blocks that appear to be off target.

     redelete:
     vec<Bool> blocks_to_delete( blocks.size( ), False );
     for ( int m = 0; m < blocks.isize( ); m++ )
     {    const double e_max = 0.1;
          if ( block_defined[0][m] && block_err_rate[0][m] >= e_max )
               blocks_to_delete[m] = True;    }
     for ( int m = 0; m < blocks.isize( ) - 1; m++ )
     {    if ( block_defined[0][m] && block_defined[0][m+1] )
          {    int over = block_stop[0][m] - block_start[0][m+1];
               if ( over <= K ) continue; // note overlaps < K are ignored

               // Does it appear that one block has a much higher error rate
               // than the other block?

               double e1 = block_err_rate[0][m], e2 = block_err_rate[0][m+1];
               const double e_min = 0.02;
               const double e_mult = 4.0;
               if ( e2 >= e_mult * Max( e1, e_min ) ) 
                    blocks_to_delete[m+1] = True;
               if ( e1 >= e_mult * Max( e2, e_min ) ) 
               {    blocks_to_delete[m] = True;    }    }    }
     if ( Sum(blocks_to_delete) > 0 )
     {    for ( int id = 0; id < N; id++ )
          {    EraseIf( block_start[id], blocks_to_delete );
               EraseIf( block_stop[id], blocks_to_delete );
               EraseIf( block_start_orig[id], blocks_to_delete );
               EraseIf( block_stop_orig[id], blocks_to_delete );
               EraseIf( block_err_rate[id], blocks_to_delete );
               EraseIf( block_defined[id], blocks_to_delete );    }
          EraseIf( blocks, blocks_to_delete );
          EraseIf( blocksx, blocks_to_delete );
          EraseIf( bpos, blocks_to_delete );
          EraseIf( bposh, blocks_to_delete );
          goto redelete;    }

     // Kill reads that have too many errors.  Currently this is done using a
     // hard threshold, which can't be right.

     {    const double max_block_err_rate = 0.06;
          for ( int id = 1; id < N; id++ )
          {    Bool bad = False;
               for ( int m = 0; m < blocks.isize( ); m++ )
                    if ( block_err_rate[id][m] > max_block_err_rate ) bad = True;
               if (bad)
               {    for ( int m = 0; m < blocks.isize( ); m++ )
                         block_defined[id][m] = False;    }    }    }

     // Prevent illegal block overlap.

     for ( int m = 0; m < blocks.isize( ) - 1; m++ )
     {    for ( int id = 0; id < (int) reads.size( ); id++ )
          {    if ( block_defined[id][m] && block_defined[id][m+1] )
               {    ForceAssertLe( block_start[id][m+1], reads[id].isize( ) ); // XXX
                    int over = block_stop[id][m] - block_start[id][m+1];
                    if ( over <= 0 ) continue;
                         
                    // Deal with weird failures.

                    // ForceAssertLt( over, blocksx[m].isize( ) );
                    if ( !( over < blocksx[m].isize( ) ) )
                    {    continue;    }

                    // Proceed with normal case.

                    blocks[m].resize( blocks[m].isize( ) - over );
                    blocksx[m].resize( blocksx[m].isize( ) - over );
                    for ( int id = 0; id < (int) reads.size( ); id++ )
                         block_stop[id][m] -= over;    }    }    }

     // Predict the position of each block on each read.  Define threaded blocks.

     vec< vec<basevector> > threads( reads.size( ) );
     vec< vec<Bool> > thread_defined( reads.size( ) );
     vec<Bool> alive( reads.size( ), False );
     vec<ho_interval> thread_range( reads.size( ) );
     if ( blocks.nonempty( ) )
     {    for ( int id = 0; id < (int) reads.size( ); id++ )
          {    threads[id].resize( blocks.size( ) - 1 );
               thread_defined[id].resize( blocks.size( ) - 1 , False );
               for ( int m = 0; m < blocks.isize( ) - 1; m++ )
               {    if ( block_defined[id][m] && block_defined[id][m+1] )
                    {    
                         // Work around weird failure.

                         if ( !( block_stop[id][m] <= block_start[id][m+1] ) )
                         {    break;    }

                         // Proceed with normal case.

                         alive[id] = True;
                         threads[id][m] = basevector( reads[id], block_stop[id][m],
                              block_start[id][m+1] - block_stop[id][m] );    
                         thread_defined[id][m] = True;    }    }
               if ( alive[id] )
               {    int thread_start = -1, thread_stop = -1;
                    for ( int j = 0; j < blocks.isize( ); j++ )
                    {    if ( thread_defined[id][j] )
                         {    thread_start = j;
                              break;    }    }
                    for ( int j = blocks.isize( ) - 2; j >= 0; j-- )
                    {    if ( thread_defined[id][j] )
                         {    thread_stop = j+1;
                              break;    }    }
                    thread_range[id] 
                         = ho_interval( thread_start, thread_stop );    }    }    }
     tb = threaded_blocks( blocks, threads, alive, thread_range );    }

int main( )
{    RunTime( );
     cout << Date( ) << ": begin" << endl;

     // Heuristics.

     const int K1 = 16;
     // ****************************************************************************
     const int K2 = 8;
     const int min_spread = 50;
     const int max_offset_diff = 100;
     const int min_start = 50;
     const int min_predicted_overlap = 1000;

     SupportedHyperBasevector shb;

     vecbasevector magic;
     vec<look_align> aligns;
     vec< vec<int> > aligns_index;


     // E. coli

     /*
     LoadLookAligns( "/wga/scr4/macro/Jan30.2014.3/all.aligns", aligns,
          aligns_index, magic.size( ) );
     magic.ReadAll( "/wga/scr4/macro/Jan30.2014.3/all.fastb" );
     BinaryReader::readFile( "/wga/scr4/macro/scs.shbv", &shb );
     */

     // Scardovia

     magic.ReadAll( "/wga/scr4/macro/scardo.Feb1/all.fastb" );
     LoadLookAligns( "/wga/scr4/macro/scardo.Feb1/all.aligns", aligns,
          aligns_index, magic.size( ) );


     int N = magic.size( );

     vecbasevector magicrc(magic);
     for ( int i = 0; i < (int) magicrc.size( ); i++ )
          magicrc[i].ReverseComplement( );
     magic.Append(magicrc);

     cout << Date( ) << ": building kmers" << endl;
     vec< triple<kmer<K1>,int,int> > kmers_plus;
     MakeKmerLookup0( magic, kmers_plus );

     /*
     cout << Date( ) << ": building kmers2" << endl;
     vec< triple<kmer<K2>,int,int> > kmers_plus2;
     MakeKmerLookup0( magic, kmers_plus2 );
     */

     cout << Date( ) << ": building matches" << endl;
     vec< triple<int,int,int> > match;
     vec< kmer<K1> > km;
     vec<int> kmult;
     vec<int> pos1s;
     for ( int i = 0; i < kmers_plus.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < kmers_plus.isize( ); j++ )
               if ( kmers_plus[j].first != kmers_plus[i].first ) break;
          for ( int k1 = i; k1 < j; k1++ )
          for ( int k2 = i; k2 < j; k2++ )
          {    if ( k2 == k1 ) continue;
               int id1 = kmers_plus[k1].second, id2 = kmers_plus[k2].second;
               if ( id1 >= N ) continue;
               int pos1 = kmers_plus[k1].third, pos2 = kmers_plus[k2].third;
               int offset = pos1 - pos2;
               match.push( id1, id2, offset );
               km.push_back( kmers_plus[i].first );
               kmult.push_back( j - i );
               pos1s.push_back(pos1);    }
          i = j - 1;    }
     // Sort(match);
     SortSync( match, km, kmult, pos1s );



     /*
     cout << Date( ) << ": building matches2" << endl;
     vec< triple<int,int,int> > match2;
     for ( int i = 0; i < kmers_plus2.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < kmers_plus2.isize( ); j++ )
               if ( kmers_plus2[j].first != kmers_plus2[i].first ) break;
          for ( int k1 = i; k1 < j; k1++ )
          for ( int k2 = i; k2 < j; k2++ )
          {    if ( k2 == k1 ) continue;
               int id1 = kmers_plus2[k1].second, id2 = kmers_plus2[k2].second;
               if ( id1 >= N ) continue;
               int pos1 = kmers_plus2[k1].third, pos2 = kmers_plus2[k2].third;
               int offset = pos1 - pos2;
               match2.push( id1, id2, offset );    }
          i = j - 1;    }
     Sort(match2);
     vec<int> start2( N + 1, -1 );
     start2[N] = match2.size( );
     for ( int i = match2.isize( ) - 1; i >= 0; i-- )
          start2[ match2[i].first ] = i;
     for ( int i = N - 1; i >= 0; i-- )
          if ( start2[i] < 0 ) start2[i] = start2[i+1];
     */



     cout << Date( ) << ": start main loop" << endl;
     vec< vec< pair<int,int> > > over(N);
     int overlaps = 0, goods = 0, bads = 0;
     for ( int i = 0; i < match.isize( ); i++ )
     {    int id1 = match[i].first, id2 = match[i].second, j;
          for ( j = i + 1; j < match.isize( ); j++ )
               if ( match[j].first != id1 || match[j].second != id2 ) break;

          int spread = 0;
          for ( int k1 = i; k1 < j; k1++ )
          for ( int k2 = i; k2 < j; k2++ )
          {    if ( Abs( match[k1].third - match[k2].third ) > max_offset_diff ) 
                    continue;
               if ( pos1s[k1] < min_start || pos1s[k2] < min_start ) continue;
               spread = Max( spread, Abs( pos1s[k1] - pos1s[k2] ) );    }

          if ( spread < min_spread ) 
          {    i = j - 1;
               continue;    }

          int offset = match[ i + (j-i)/2 ].third;
          over[id1].push( id2, offset );

          overlaps++;

          cout << "\n";

          /*
          int support = 0;
          for ( int l = start2[id1]; l < start2[id1+1]; l++ )
               if ( match2[l].first == id1 && match2[l].second == id2 ) support++;
          PRINT(support);
          */

          /*
          int id1x(id1), id2x(id2);
          if ( id1 >= N ) id1x -= N;
          if ( id2 >= N ) id2x -= N;
          if ( aligns_index[id1x].solo( ) && aligns_index[id2x].solo( ) )
          {    const look_align& la1 = aligns[ aligns_index[id1x][0] ];
               const look_align& la2 = aligns[ aligns_index[id2x][0] ];
               int e1 = la1.target_id, e2 = la2.target_id;
               if ( la1.pos2( ) > 0 && la2.pos2( ) > 0
                    && la1.Pos2( ) < shb.EdgeLengthBases(e1)
                    && la2.Pos2( ) < shb.EdgeLengthBases(e2) )
               {    Bool rc1 = ( id1 >= N ), rc2 = ( id2 >= N );
                    if (rc1) e1 = shb.Inv(e1);
                    if (rc2) e2 = shb.Inv(e2);
                    int start1 = la1.pos2( ), stop1 = la1.Pos2( );
                    int start2 = la2.pos2( ), stop2 = la2.Pos2( );
                    if (rc1)
                    {    int start1x = shb.EdgeLengthBases(e1) - stop1;
                         int stop1x = shb.EdgeLengthBases(e1) - start1;
                         start1 = start1x;
                         stop1 = stop1x;    }
                    if (rc2)
                    {    int start2x = shb.EdgeLengthBases(e2) - stop2;
                         int stop2x = shb.EdgeLengthBases(e2) - start2;
                         start2 = start2x;
                         stop2 = stop2x;    }
                    Bool bad = False;
                    if ( e1 != e2 )
                    {    cout << "different targets" << endl;
                         bad = True;    }
                    else 
                    {    int overlap = IntervalOverlap(start1, stop1, start2, stop2);
                         if ( overlap <= 0 ) 
                         {    cout << "nothing" << endl;
                              bad = True;    }
                         else 
                         {    cout << "overlap" << endl;
                              goods++;     }    }
                    if (bad) 
                    {    bads++;    
                         for ( int k = i; k < j; k++ )
                         {    String the_kmer = km[k].ToString( );
                              int mult = kmult[k];
                              int pos1 = pos1s[k];
                              PRINT3( the_kmer, mult, pos1 );    }    }    }    }
          */
               
          for ( int k = i; k < j; k++ )
          {    int l;
               for ( l = k; l < j; l++ )
                    if ( match[l].third != match[k].third ) break;
               int count = l - k;
               int offset = match[k].third;
               cout << "id1 = ";
               if ( id1 < N ) cout << id1;
               else cout << id1 - N << "'";
               cout << ", id2 = ";
               if ( id2 < N ) cout << id2;
               else cout << id2 - N << "'";
               cout << ", ";
               PRINT3( offset, count, spread );

               /*
               for ( l = k; l < j; l++ )
               {    if ( match[l].third != match[k].third ) break;
                    String the_kmer = km[l].ToString( );
                    int mult = kmult[l];
                    int pos1 = pos1s[l];
                    PRINT3( the_kmer, mult, pos1 );    }
               */

               k = l - 1;    }
          i = j - 1;    }

     cout << "\n";
     PRINT3( overlaps, goods, bads );
     cout << "bads/goods = " << PERCENT_RATIO( 3, bads, goods ) << endl;
     cout << Date( ) << ": done with main loop" << endl;

     // Generate transitive overlaps.

     vec<String> reports(N);
     #pragma omp parallel for
     for ( int id1 = 0; id1 < N; id1++ )
     {    
          if ( id1 != 1006 ) continue; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

          vec< pair<int,int> > all_over = over[id1];
          ostringstream out;
          out << "\noverlaps with id1 = " << id1 
               << ", len = " << magic[id1].size( ) << "\n";
          out << "direct:\n";
          vec<int> ids;
          ids.push_back(id1);
          for ( int i = 0; i < over[id1].isize( ); i++ )
          {    int id2 = over[id1][i].first;
               ids.push_back(id2);
               out << "id2 = ";
               if ( id2 < N ) out << id2;
               else out << id2 - N << "'";
               out << ", len = " << magic[id2].size( )
                    << ", offset = " << over[id1][i].second << "\n";    }
          Sort(ids);
          vec<int> ids_plus = ids;
          out << "transitive:\n";
          vec< pair<int,int> > trans;
          for ( int i = 0; i < over[id1].isize( ); i++ )
          {    int id2 = over[id1][i].first;
               int offset2 = over[id1][i].second;
               Bool rc2 = False;
               if ( id2 >= N )
               {    id2 -= N;
                    rc2 = True;    }
               for ( int j = 0; j < over[id2].isize( ); j++ )
               {    int id3 = over[id2][j].first;
                    int offset3 = over[id2][j].second;
                    if (rc2)
                    {    if ( id3 < N ) id3 += N;
                         else id3 -= N;    
                         offset3 = -offset3 
                              - magic[id3].isize( ) + magic[id2].isize( );    }
                    trans.push( id3, offset2 + offset3 );    }    }
          Sort(trans);
          for ( int t = 0; t < trans.isize( ); t++ )
          {    int u;
               for ( u = t + 1; u < trans.isize( ); u++ )
                    if ( trans[u].first != trans[t].first ) break;
               int offset = 0;
               for ( int l = t; l < u; l++ )
                    offset += trans[l].second;
               offset /= (u-t);
               int id3 = trans[t].first;
               if ( !BinMember( ids, id3 ) )
               {    ids_plus.push_back(id3);
                    int pred = IntervalOverlap( 0, magic[id1].isize( ),
                         offset, offset + magic[id3].isize( ) );
                    if ( pred < min_predicted_overlap ) continue;
                    out << "id3 = ";
                    if ( id3 < N ) out << id3;
                    else out << id3 - N << "'";
                    out << ", len = " << magic[id3].size( )
                         << ", offset = " << offset 
                         << ", predicted overlap = " << pred << "\n";
                    all_over.push( id3, pred );    }
               t = u - 1;    }

          // Make gang and alignments.

          vecbasevector gang;
          for ( int i = 0; i < all_over.isize( ); i++ )
               gang.push_back( magic[ all_over[i].first ] );





     // Generate multiple alignment.  Note redundant generation of pairwise alignments.

     #pragma omp critical
     {    cout << "\nmultiple alignment\n\n";
          const double del_rate = 0.05;
          const double ins_rate = 0.02;
          const double sub_rate = 0.05;
          Scorer scorer( sub_rate, del_rate, ins_rate );
          VecUCharVec multi;
          vecbasevector gang;
          gang.push_back( magic[id1] );
          vec<align> aligns;
          aligns.push_back( align( ) );
          const int bandwidth = 200;
          for ( int i = 0; i < all_over.isize( ); i++ )
          {    gang.push_back( magic[ all_over[i].first ] );
               align x;
               int errors;
               SmithWatBandedA( magic[ all_over[i].first ],  magic[id1] ,
                    -all_over[i].second, bandwidth, x, errors, 0, 1, 1 );    
               x.Flip( );
               aligns.push_back(x);    }
          AlignFriendsToFounder( gang, 0, aligns, scorer, &multi );

          // Delete columns having only one base, and at least five entries.

          vec<Bool> to_delete( multi[0].size( ), False );
          for ( int j = 0; j < (int) multi[0].size( ); j++ )
          {    int non45 = 0;
               for ( int i = 0; i < (int) multi.size( ); i++ )
                    if ( multi[i][j] != 4 && multi[i][j] != 5 ) non45++;
               if ( non45 == 1 && multi.size( ) >= 5 ) to_delete[j] = True;    }
          for ( int i = 0; i < (int) multi.size( ); i++ )
          {    SerfVec<uchar> m = multi[i];
               vec<uchar> mx;
               for ( int l = 0; l < (int) m.size( ); l++ )
                    mx.push_back( m[l] );
               EraseIf( mx, to_delete );
               multi[i].resize(0);
               for ( int l = 0; l < mx.isize( ); l++ )
                    multi[i].push_back( mx[l] );    }

          // Print the multi-alignment.

          const int width = 80;
          for ( int js = 0; js < (int) multi[0].size( ); js += width )
          {    for ( int i = 0; i < (int) multi.size( ); i++ )
               {    Bool nothing = True;
                    for ( int j = js; j < Min( (int) multi[0].size( ), js + width ); j++ )
                    {    uchar c = multi[i][j];
                         if ( c != 5 ) nothing = False;    }
                    if ( !nothing )
                    {    for ( int j = js; j < Min( (int) multi[0].size( ), js + width ); 
                              j++ )
                         {    uchar c = multi[i][j];
                              if ( c < 4 ) cout << as_base(c);
                              if ( c == 4 ) cout << "-";
                              if ( c == 5 ) cout << "=";    }
                         cout << "\n";    }    }
               cout << "\n";    }    }





          vec< pair<int,char> > subs;
          for ( int i = 0; i < all_over.isize( ); i++ )
          {    
               align x;
               int errors;
               SmithWatBandedA( magic[ all_over[i].first ], magic[id1],
                    -all_over[i].second, 200, x, errors, 0, 1, 1 );    
               if ( id1 == 1006 )
               {
               out << "\nalignment to read " << all_over[i].first << endl;
               PrintVisualAlignment( True, out,
                    magic[ all_over[i].first ], magic[id1], x );
               vec<ho_interval> perf1;
               x.PerfectIntervals1( magic[ all_over[i].first ], magic[id1], perf1 );
               int M = 0;
               for ( int l = 0; l < perf1.isize( ); l++ )
                    M = Max( M, perf1[l].Length( ) );
               out << "max perfect match = " << M << "\n";

               basevector rd1 = magic[ all_over[i].first ];
               basevector rd2 = magic[id1];

               int p1 = x.pos1( ), p2 = x.pos2( );
               for ( int j = 0; j < x.Nblocks( ); j++ )
               {    if ( x.Gaps(j) > 0 ) 
                    {    out << "EDIT: delete " 
                              << x.Gaps(j) << " bases at " << p2 << endl;
                         p2 += x.Gaps(j);    }
                    if ( x.Gaps(j) < 0 ) 
                    {    out << "EDIT: insert ";
                         for ( int l = 0; l < -x.Gaps(j); l++ )
                              out << as_base( rd1[ p1 + l ] );
                         out << " at " << p2 << endl;
                         p1 -= x.Gaps(j);    }
                    for ( int l = 0; l < x.Lengths(j); l++ )
                    {    if ( rd1[p1] != rd2[p2] )
                         {    out << "EDIT: change to " << as_base( rd1[p1] )
                                   << " at " << p2 << endl;
                              subs.push( p2, rd1[p1] );    }
                         ++p1; ++p2;    }    }    }    }

          Sort(subs);
          basevector mod( magic[id1] );
          for ( int l = 0; l < subs.isize( ); l++ )
          {    int m = subs.NextDiff(l);
               if ( m - l >= 4 )
                    mod.Set( subs[l].first, subs[l].second );
               l = m - 1;    }
          cout << out.str( );
          mod.Print( cout, "mod" );
          Scram(0);
               
          reports[id1] = out.str( );
          continue;

          // Make gang.  Following paths/long/ultra/GetFriendsAndAlignsInitial.cc.

          vec< vec< pair<int,int> > > a;

          vec< triple<kmer<K2>,int,int> > kmers_plus;
          MakeKmerLookup0SingleX( gang, kmers_plus );
          int N = gang.size( ); // argh, overusing this symbol
          vec< vec< pair<int,int> > > offsets(N);
          for ( int64_t i = 0; i < (int64_t) kmers_plus.size( ); i++ )
          {    int64_t j, z1, z2;
               for ( j = i + 1; j < (int64_t) kmers_plus.size( ); j++ )
                    if ( kmers_plus[j].first != kmers_plus[i].first ) break;
               for ( z1 = i; z1 < j; z1++ )
                    if ( kmers_plus[z1].second == 0 ) break;
               for ( z2 = z1; z2 < j; z2++ )
                    if ( kmers_plus[z2].second != 0 ) break;
               for ( int64_t z = z1; z < z2; z++ )
               {    int pos1 = kmers_plus[z].third;
                    for ( int64_t k = i; k < j; k++ )
                    {    int id2 = kmers_plus[k].second;
                         if ( id2 == 0 ) continue;
                         int pos2 = kmers_plus[k].third; // position on read id2
                              offsets[id2].push( pos1-pos2, pos1 );    }    }
               i = j - 1;    }
          for ( int id = 1; id < N; id++ )
               Sort( offsets[id] );

          vec<align> aligns( gang.size( ) );
          vec<Bool> accepted( N, False );
          vec<int> errs(N, -1);
          vec<ho_interval> ext1(N);
          {    int mh = 0;
               for ( int id = 0; id < N; id++ )
                    if ( id != 0 ) mh = Max( mh, offsets[id].isize( ) );
               vec<int> pos1_count( N );
               for ( int id = 0; id < N; id++ )
               {    if ( id == 0 ) continue;
                    vec<int> pos1;
                    for ( int j = 0; j < offsets[id].isize( ); j++ )
                         pos1.push_back( offsets[id][j].second );
                    UniqueSort(pos1);
                    pos1_count[id] = pos1.size( );    }
               int max_pos1_count = 0;
               for ( int id = 0; id < N; id++ )
               {    if ( id != 0 ) 
                         max_pos1_count = Max( max_pos1_count, pos1_count[id] );    }


               a.clear_and_resize(N);
               for ( int id = 0; id < N; id++ )
               {    if ( id == 0 ) continue;
                    if ( pos1_count[id] >= max_pos1_count/2 
                         && offsets[id].nonempty( ) )
                    {    accepted[id] = True;

                         // Form offsets into groups, breaking whenever there is a
                         // > max_sep = 20 separation in offsets.

                         vec<int> ostarts;
                         ostarts.push_back(0);
                         const int max_sep = 20;
                         for ( int j = 0; j < offsets[id].isize( ) - 1; j++ )
                         {    if ( offsets[id][j+1].first 
                                   > offsets[id][j].first + 20 )
                              {    ostarts.push_back(j+1);    }    }
                         ostarts.push_back( offsets[id].size( ) );

                         // Compute the size of each group, as measured by its number
                         // of distinct rpos2 values.  Find the largest group.

                         int gp_max = 0, gp_best = -1;
                         for ( int j = 0; j < ostarts.isize( ) - 1; j++ )
                         {    vec<int> rpos2;
                              for ( int l = ostarts[j]; l < ostarts[j+1]; l++ )
                                   rpos2.push_back( offsets[id][l].second );
                              UniqueSort(rpos2);
                              if ( rpos2.isize( ) > gp_max )
                              {    gp_max = rpos2.size( );
                                   gp_best = j;    }    }

                         // Set offset to the median within the winning group, and
                         // set bandwidth to span the group, but no more than 
                         // max_bandwidth = 200;

                         int offset, bandwidth;
                         {    int start = ostarts[gp_best]; 
                              int stop = ostarts[gp_best+1];
                              int mid = start + (stop-start)/2;
                              offset = offsets[id][mid].first;
                              const int bw_add = 12;
                              bandwidth = Max( offset - offsets[id][start].first,
                                   offsets[id][stop-1].first - offset ) + bw_add;
                              const int max_bandwidth = 200;
                              bandwidth = Min( bandwidth, max_bandwidth );    }
     
                         if ( bandwidth < 0 ) // SHOULD NOT HAPPEN!!!!!!!!!!!!!!!!!!
                              bandwidth = 0;
     
                         // Align.

                         align x;
                         int errors;
                         if ( id1 == 1006 ) bandwidth += 500;
                         SmithWatBandedA( gang[0], gang[id], offset, bandwidth, x,
                              errors, 0, 1, 1 );    

                         if ( id1 == 1006 )
                         {
                         out << "\nalignment to read " << id << endl;
                         PrintVisualAlignment( True, out, gang[0], gang[id], x );
                         vec<ho_interval> perf1;
                         x.PerfectIntervals1( gang[0], gang[id], perf1 );
                         int M = 0;
                         for ( int l = 0; l < perf1.isize( ); l++ )
                              M = Max( M, perf1[l].Length( ) );
                         out << "max perfect match = " << M << "\n";
                         }

                         errs[id] = errors;
                         ext1[id] = x.Extent1( );
                         aligns[id] = x;
                         vec<ho_interval> p1, p2;
                         x.PerfectIntervals1( gang[0], gang[id], p1 );
                         x.PerfectIntervals2( gang[0], gang[id], p2 );    
                         for ( int j = 0; j < p1.isize( ); j++ )
                         {    const ho_interval &h1 = p1[j], &h2 = p2[j];
                              for ( int l = h1.Start( ); l <= h1.Stop( ) - K2; l++ )
                         {    a[id].push( l + h2.Start( ) - h1.Start( ), 
                                   l );    }    }    }    }    }

          for ( int p = 0; p <= gang[0].isize( ) - K2; p++ )
               a[0].push( p, p );

          // Filter out reads based on error rate.  Note that this doesn't actually
          // remove the reads yet (which might be better).

          vec<Bool> evil( gang.size( ), False );
          {    vec<Bool> e_to_delete( N, False );
               vec<int> eids( N, vec<int>::IDENTITY );
               for ( int i = 0; i < N; i++ )
                    if ( errs[i] < 0 ) e_to_delete[i] = True;
               EraseIf( errs, e_to_delete );
               EraseIf( ext1, e_to_delete );
               EraseIf( eids, e_to_delete );
               vec<double> erate( errs.size( ) );
               for ( int j = 0; j < errs.isize( ); j++ )
                    erate[j] = double(errs[j]) / double( ext1[j].Length( ) );
               SortSync( erate, errs, ext1, eids );
     
               // Walk through errs until the founder is covered.

               vec<ho_interval> cov;
               int ep;
               for ( ep = 0; ep < ext1.isize( ); ep++ )
               {    cov.push_back( ext1[ep] );
                    if ( TotalCovered(cov) == gang[0].isize( ) ) break;    }
               if ( ep < ext1.isize( ) )
               {    
                    // Filter out reads that have too high an error rate.
     
                    const double dev_mult = 3.0;
                    double add = dev_mult * erate[ep] * sqrt(errs[ep]) / errs[ep];
                    for ( int j = 0; j < errs.isize( ); j++ )
                    {    if ( erate[j] > erate[ep] + add )
                              evil[ eids[j] ] = True;    }    }    }

          // Filter out reads that appear not to belong.

          vec<double> evil_score( gang.size( ), 1 );
          VecUCharVec multi;
          vec<Bool> funkyc( gang[0].size( ), False );
          {
               // Create substitution-only multiple alignment.

               unsigned char EMPTY_CELL = 5;
               int nrows = gang.size( ), ncols = gang[0].size();
               {
                   multi.resize(nrows);
                   for ( int id = 0; id < nrows; id++ )
                   {    multi[id].resize( ncols, EMPTY_CELL );
                        if ( id == 0 )
                        {    for ( int j = 0; j < ncols; j++ )
                                  multi[id][j] = gang[0][j];    }
                        else
                        {    const align& a = aligns[id];
                             int p1 = a.pos1( ), p2 = a.pos2( );
                             for ( int j = 0; j < a.Nblocks( ); j++ )
                             {    if ( a.Gaps(j) > 0 ) p2 += a.Gaps(j);
                                  if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
                                  for ( int x = 0; x < a.Lengths(j); x++ )
                                  {    multi[id][p1] = gang[id][p2];
                                       ++p1; ++p2;    }    }    }    }    }

               // Look for funky columns.
     
               // For haploid genome, we may want the threshold to be strict so that
               // to reduce false negatives. However it is necessary to increase the
               // threshold for polymorphic genome.
               // double eps = 1.0e-40;

               // ******************************************************************
               /*
               const double eps = 1.0e-20;
               VecUCharVec F(nrows);
               for ( int c = 0; c < ncols; c++ )
               {    Bool funky = True;
                    vec<int> x(5, 0); // (A,C,G,T,-)
                    for ( int r = 0; r < nrows; r++ )
                         if ( multi[r][c] != EMPTY_CELL ) x[ multi[r][c] ]++;
                    if ( logc.verb[ "FRIEND" ] >= 2 )
                    {    out << "\n";
                         PRINT6_TO( out, c, x[0], x[1], x[2], x[3], x[4] );    }
                    int N = Sum(x);
                    double worst = 0.0;
                    for ( int j = 0; j < 5; j++ )
                    {    Bool funkyj = False;
                         double w = 1.0;
                         for ( int k = 0; k < 5; k++ )
                         {    int n = x[k];
                              if ( k == j || n == 0 ) continue;
                              double p;
                              if ( j < 4 )
                              {    if ( k < 4 ) p = error_model.GetSubRate( )/3.0;
                                   else p = error_model.GetDelRate( );    }
                              else p = error_model.GetInsRate( )/3.0;
                              double q = BinomialProbCum( 1.0-p, N, N-n );
                              if ( q < w ) w = q;
                              if ( q <= eps )
                              {    funkyj = True;
                                   // break;    
                                   }    }
                         if ( w > worst ) worst = w;
                         if ( !funkyj )
                         {    funky = False;
                              break;    }    }
                    if ( !funky ) continue;
                    funkyc[c] = True;
                    if ( logc.verb[ "FRIEND" ] >= 1 )
                    {    out << "column " << c << " looks funky, score = " 
                              << worst << "\n";    }
                    for ( int r = 0; r < nrows; r++ )
                         F[r].push_back( multi[r][c] );    }

               // Look for reads that appear not to belong.

               const double pd = 0.05;
               const double eps2 = 0.000001; // 10^-6
               for ( int j = 0; j < nrows; j++ )
               {    if ( j == 0 ) continue;
                    if ( !accepted[j] ) continue;
                    int total = 0, diffs = 0;
                    for ( size_t c = 0; c < F[j].size( ); c++ )
                    {    if ( F[j][c] == EMPTY_CELL ) continue;
                         total++;
                         if ( F[j][c] != F[0][c] ) diffs++;    }
                    if ( diffs == 0 ) continue;
                    double q = 1.0 - BinomialSum( total, diffs-1, pd );
                    evil_score[j] = q;
                    if ( q <= eps2 )
                    {    evil[j] = True;    }    }
               */
                    }

          // Delete reads for which number of alignments is less than half the max.
          // Print gang.

          int M = 0;
          for ( int id = 1; id < N; id++ )
               M = Max( M, a[id].isize( ) );
          vec<Bool> to_delete( N, False );
          vec< vec<String> > gang_rows;
          vec< pair<int,int> > gang_start;
          vec<int> colsp;

          int count = 0, true_friends = 0, false_friends = 0;
          for ( int id = 1; id < N; id++ )
          {    if ( a[id].isize( ) < M/2 || a[id].empty( ) /* || evil[id] */ ) 
                    to_delete[id] = True;
               if ( evil[id] ) to_delete[id] = True;    }
          gang.EraseIf(to_delete);
          out << "final gang size = " << gang.size( ) << endl;
          EraseIf( a, to_delete );
          // EraseIf( friend_ids, to_delete );
          // if ( p_friends != 0 )  *p_friends = friend_ids;

          // Make threaded blocks, then build "corrected read", and trim K
          // bases off its ends.  Now following ultra/MakeBlocks.cc.

          threaded_blocks tb;
          MakeBlocksX<K2>( gang, a, tb, out );
          PRINT3_TO( out, tb.NBlocks( ), tb.NGaps( ), tb.NReads( ) );
          for ( int l = 0; l < tb.NBlocks( ); l++ )
               PRINT3_TO( out, l, tb.Block(l).size( ), tb.Block(l).ToString( ) );

          // ***********************************************************************
          ConsensusScoreModel error_model( 0.05, 0.05, 0.05 );


          long_logging logc( "", "" );
          logc.STATUS_LOGGING = False;
          logc.MIN_LOGGING = False;
          ref_data ref;
          vec<ref_loc> readlocs;
          long_logging_control log_control( ref, &readlocs, "", "" );
          long_heuristics heur( "" );

          efasta r = tb.MakeCorrectedRead(
               error_model, out, heur, log_control, logc );

          out << "corrected read = \n";
          r.Print(out);

          reports[id1] = out.str( );

          }

     for ( int i = 0; i < reports.isize( ); i++ )
          cout << reports[i];

     cout << "\n" << Date( ) << ": done" << endl;    }
