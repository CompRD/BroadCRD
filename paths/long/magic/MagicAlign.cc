///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Before running this, should have run from ~jaffe/crd:
// Porkula N=-1 OUT_FASTB= ~/crd/woof.fastb VISUAL=True >&  ~/crd/Porkula.out

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "Equiv.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "PrintAlignment.h"
#include "math/Functions.h"
#include "pairwise_aligners/ClusterAligner.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/CreateGenome.h"

int main( )
{
     RunTime( );

     // Experimental constants.

     int NNN = 2000; // temporary cap on number of reads
     int KEY = 133; // key read

     // Heuristics.

     const int MAX_OFFSET_DIFF = 1000;
     const int MIN_CLUSTER_SIZE = 20;
     const int MAX_FREQ = 100;
     const int BW_ADD = 300;
     const int MISMATCH_PENALTY = 2;
     const int GAP_PENALTY = 3;
     const int MAX_FUDGE = 10;
     const int K = 12;
          
     // Load reads.

     cout << Date( ) << ": loading reads" << endl;
     vecbasevector all( "/wga/dev/jaffe/BroadCRD/woof.fastb" );
     all.resize(NNN); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     int nreads = all.size( );
     PRINT(nreads);
     vecbasevector all_rc(all);
     for ( int i = 0; i < (int) all.size( ); i++ )
          all_rc[i].ReverseComplement( );
     all.Append(all_rc);

     // Fetch read locations on the genome.

     vec<int> loc_g(nreads, -1), loc_start(nreads), loc_stop(nreads);
     fast_ifstream in( "/wga/dev/jaffe/BroadCRD/Porkula.out" );
     String line;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          if ( line.Contains( "LOC", 0 ) )
          {    int id = line.Between( "LOC: ", " " ).Int( );
               if ( id >= nreads ) continue; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               int g = line.Between( "at ", "." ).Int( );
               int start = line.Between( ".", "-" ).Int( );
               int stop = line.After( "-" ).Int( );
               loc_g[id] = g;
               loc_start[id] = start;
               loc_stop[id] = stop;    }    }

     // Make lookup table for genome.  Compress homopolymers.

     cout << Date( ) << ": making lookup table for genome" << endl;
     vecbasevector GG( "/wga/dev/references/Escherichia_coli/genome.fastb" );
     for ( int g = 0; g < (int) GG.size( ); g++ )
     {    basevector x;
          for ( int j = 0; j < GG[g].isize( ); j++ )
          {    if ( j >= 5 && GG[g][j] == GG[g][j-1] && GG[g][j] == GG[g][j-2]
                    && GG[g][j] == GG[g][j-3] && GG[g][j] == GG[g][j-4]
                    && GG[g][j] == GG[g][j-5] )
               {    continue;     }
               x.push_back( GG[g][j] );    }
          GG[g] = x;    }

     // Print test interval.

     cout << "\n";
     basevector Z;
     Z.SetToSubOf( GG[0], 72000, 1000 );
     Z.Print( cout, "Z" );
     cout << "\n";

     vecbasevector genome(GG);
     vecbasevector genome_rc(genome);
     for ( int g = 0; g < (int) genome.size( ); g++ )
          genome_rc[g].ReverseComplement( );
     genome.Append(genome_rc);
     vec< triple<kmer<K>,int,int> > gkmers_plus;
     MakeKmerLookup0( genome, gkmers_plus );

     // Make kmer lookup table.

     cout << Date( ) << ": making kmers" << endl;
     vec< triple<kmer<K>,int,int> > kmers_plus;
     MakeKmerLookup0( all, kmers_plus );

     // Find hits.

     cout << Date( ) << ": finding hits" << endl;
     vec< triple<int,int,int> > hits;
     for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
     {    int64_t j;
          for ( j = i + 1; j < kmers_plus.jsize( ); j++ )
               if ( kmers_plus[j].first != kmers_plus[i].first ) break;
          vec< triple<int,int,int> > x;
          if ( j - i <= MAX_FREQ )
          {    for ( int64_t k1 = i; k1 < j; k1++ )
               {    int id1 = kmers_plus[k1].second;
                    if ( id1 >= nreads ) continue;
                    int pos1 = kmers_plus[k1].third;
                    for ( int64_t k2 = i; k2 < j; k2++ )
                    {    int id2 = kmers_plus[k2].second;
                         if ( id1 >= id2 ) continue;
                         int pos2 = kmers_plus[k2].third;
                         hits.push( id1, id2, pos2 - pos1 );    }    }    }
          i = j - 1;    }
     PRINT( hits.size( ) );

     // Sort hits.

     cout << Date( ) << ": sorting hits" << endl;
     ParallelSort(hits);

     // Find clusters.

     vec< quad<int,int,int,int> > clusters;
     cout << Date( ) << ": finding clusters" << endl;
     for ( int64_t i = 0; i < hits.jsize( ); i++ )
     {    int64_t j;
          int id1 = hits[i].first, id2 = hits[i].second;
          for ( j = i + 1; j < hits.jsize( ); j++ )
          {    if ( hits[j].first != id1 ) break;
               if ( hits[j].second != id2 ) break;
               if ( hits[j].third - hits[i].third > MAX_OFFSET_DIFF ) break;    }
          if ( j - i >= MIN_CLUSTER_SIZE ) 
          {    int id1 = hits[i].first, id2 = hits[i].second;
               int offset_low = hits[i].third;
               int offset_high = hits[j-1].third;
               int offset = (offset_low + offset_high) / 2;
               int bandwidth = (offset_high - offset_low) / 2;
               clusters.push( id1, id2, offset, bandwidth );    }
          i = j - 1;    }
     PRINT( clusters.size( ) );

     // Align.

     cout << Date( ) << ": aligning" << endl;
     vec< triple<int,int,align> > aligns;
     #pragma omp parallel for
     for ( int64_t i = 0; i < clusters.jsize( ); i++ )
     {    int id1 = clusters[i].first, id2 = clusters[i].second;
          int offset = clusters[i].third, bandwidth = clusters[i].fourth;
          int errsx;
          align a;
          SmithWatBandedA( all[id1], all[id2], -offset, bandwidth + BW_ADD, 
               a, errsx, 0, MISMATCH_PENALTY, GAP_PENALTY );
          vec<ho_interval> p1;
          a.PerfectIntervals1( all[id1], all[id2], p1 );
          int ngood = 0;
          for ( int l = 0; l < p1.isize( ); l++ )
               if ( p1[l].Length( ) >= K ) ngood += p1[l].Length( ) - K + 1;
          if ( ngood < MIN_CLUSTER_SIZE ) continue;
          #pragma omp critical
          {    // PRINT4( id1, id2, offset, ngood );
               // PrintVisualAlignment( True, cout, all[id1], all[id2], a );
               aligns.push( id1, id2, a );    }    }
     PRINT( aligns.size( ) );

     // Sort aligns.

     cout << Date( ) << ": sorting alignments" << endl;
     ParallelUniqueSort(aligns);

     // Summarize alignments.

     vec< vec< pair<int,align> > > alignsx(nreads);
     for ( int i = 0; i < aligns.isize( ); i++ )
          alignsx[ aligns[i].first ].push( aligns[i].second, aligns[i].third );

     // Add reverse alignments.

     alignsx.resize( 2*nreads );
     for ( int id1 = nreads; id1 < 2*nreads; id1++ )
     {    int id1r = id1 - nreads;
          for ( int i = 0; i < alignsx[id1r].isize( ); i++ )
          {    int id2r = alignsx[id1r][i].first;
               align a = alignsx[id1r][i].second;
               int id2 = id2r;
               if ( id2 >= nreads ) id2 -= nreads;
               else id2 += nreads;
               a.ReverseThis( all[id1].size( ), all[id2].size( ) );
               alignsx[id1].push( id2, a );    }    }

     for ( int id1 = 0; id1 < nreads; id1++ )
     {    cout << "read " << id1 << " aligned to";
          for ( int j = 0; j < alignsx[id1].isize( ); j++ )
          {    int id2 = alignsx[id1][j].first;
               cout << " " << id2;    }
          cout << "\n";

          // Build graph for key read.

          if ( id1 == KEY )
          {    vec<int> ids = {id1};
               for ( int j = 0; j < alignsx[id1].isize( ); j++ )
               {    int id2 = alignsx[id1][j].first;
                    ids.push_back(id2);    }
               UniqueSort(ids);
               vec<String> kmers;
               int N = 0; // total kmers
               int K = 5;
               for ( int i = 0; i < ids.isize( ); i++ )
               {    N += all[ ids[i] ].isize( ) - K + 1;
                    for ( int j = 0; j < all[ ids[i] ].isize( ) - K + 1; j++ )
                    {    basevector b;
                         b.SetToSubOf( all[ ids[i] ], j, K );
                         kmers.push_back( b.ToString( ) );    }    }

               vec<vec<int>> locs( kmers.size( ) );
               const int L = 12;
               VecIntPairVec Glocs;
               CreateGlocs( genome, L, Glocs );
               for ( int i = 0; i < ids.isize( ); i++ )
               {    vec<look_align> aligns;
                    int BW_ADD = 300;
                    const int MIN_CLUSTER = 20;
                    const int MAX_OFFSET_DIFF = 1000;
                    const int MISMATCH_PENALTY = 2;
                    const int GAP_PENALTY = 3;
                    Bool FW_ONLY = True;
                    const basevector& q = all[ ids[i] ];
                    ClusterAligner( q, genome, L, Glocs, aligns, FW_ONLY, 
                         BW_ADD, MIN_CLUSTER, MAX_OFFSET_DIFF, 
                         MISMATCH_PENALTY, GAP_PENALTY );
                    if ( aligns.empty( ) ) continue;
                    // NOTE NOT TRACKING GENOME CONTIG ID!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    int g = aligns[0].target_id;
                    const align& a = aligns[0].a;
                    int p1 = a.pos1( ), p2 = a.pos2( );
                    // PRINT3( ids[i], p1, p2 ); // YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
                    // PrintVisualAlignment( True, cout, q, genome[g], a ); // YYYYY
                    cout << "read " << ids[i] << " aligned to genome at "
                         << g << "." << a.pos2( ) << "-" << a.Pos2( ) << endl;
                    for ( int j = 0; j < a.Nblocks( ); j++ ) 
                    {    if ( a.Gaps(j) > 0 ) p2 += a.Gaps(j);
                         if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
                         for ( int x = 0; x < a.Lengths(j); x++ ) 
                         {    if ( x <= a.Lengths(j) - K )
                              {    Bool mismatch = False;
                                   for ( int i = 0; i < K; i++ )
                                   {    if ( q[p1+i] != genome[g][p2+i] )
                                             mismatch = True;    }
                                   if ( !mismatch )
                                   {    int n = p1;
                                        for ( int l = 0; l < i; l++ )
                                             n += all[ ids[l] ].isize( ) - K + 1;
                                        locs[n].push_back(p2);    }    }
                              ++p1;
                              ++p2;    }    }    }

               // String query = "CGCTGGCT";
               // String query = "TTAAGCTAAA"; // 72771
               String query = "CGCTGGCTGG"; // 72672
               for ( int i = 0; i < ids.isize( ); i++ )
               {    int idi = ids[i];
                    String s = all[ ids[i] ].ToString( );
                    for ( int j = 0; j < s.isize( ); j++ )
                    {    if ( s.Contains( query, j ) )
                         {    cout << "FOUND " << s.substr( j-5, 5 )
                                   << "|" << query << "|" 
                                   << s.substr( j+query.isize( ), 5 )
                                   << " AT " << idi << "." << j 
                                   << endl;    }    }    }

               // Build equivalence relation by joining reads along alignments.

               int pi1 = BinPosition( ids, id1 );
               equiv_rel e(N);
               for ( int i = 0; i < ids.isize( ); i++ )
               {    int idi = ids[i];

                    for ( int k = 0; k < alignsx[idi].isize( ); k++ )
                    {    int idj = alignsx[idi][k].first;
                         int j = BinPosition( ids, idj );
                         if ( j < 0 ) continue;
                         if ( !( i < j ) ) continue;
                         align a = alignsx[idi][k].second;
                         cout << "\nalignment of " << idi << " to " << idj 
                              << ", offset = " << a.pos1( ) - a.pos2( ) << endl;
                         PrintVisualAlignment( True, cout, all[idi], // YYYYYYYYYYYY
                              all[idj], a ); // YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
                         vec<ho_interval> perf1, perf2;
                         a.PerfectIntervals1( all[idi], all[idj], perf1 );
                         a.PerfectIntervals2( all[idi], all[idj], perf2 );

                         // Extend perfect intervals.

                         for ( int i = 0; i < perf1.isize( ); i++ )
                         {    while(1)
                              {    if ( perf1[i].Start( ) == 0 ) break;
                                   if ( perf2[i].Start( ) == 0 ) break;
                                   if ( all[idi][ perf1[i].Start( ) - 1 ]
                                        != all[idj][ perf2[i].Start( ) - 1 ] )
                                   {    break;    }
                                   perf1[i].AddToStart(-1);
                                   perf2[i].AddToStart(-1);    }
                              while(1)
                              {    if ( perf1[i].Stop( ) == all[idi].isize( ) )
                                        break;
                                   if ( perf2[i].Stop( ) == all[idj].isize( ) )
                                        break;
                                   if ( all[idi][ perf1[i].Stop( ) ]
                                        != all[idj][ perf2[i].Stop( ) ] )
                                   {    break;    }
                                   perf1[i].AddToStop(1);
                                   perf2[i].AddToStop(1);    }    }

                         int low, high;
                         for ( low = 0; low < perf1.isize( ); low++ )
                              if ( perf1[low].Length( ) >= 10 ) break;
                         for ( high = perf1.isize( ) - 1; high >= 0; high-- )
                              if ( perf1[high].Length( ) >= 10 ) break;
                         for ( int m = low; m <= high; m++ )
                         for ( int u = 0; u <= perf1[m].Length( ) - K; u++ )
                         {    int p1 = perf1[m].Start( ) + u;
                              int p2 = perf2[m].Start( ) + u;
                              int n1 = p1, n2 = p2;
                              for ( int l = 0; l < i; l++ )
                                   n1 += all[ids[l]].isize( ) - K + 1;
                              for ( int l = 0; l < j; l++ )
                                   n2 += all[ids[l]].isize( ) - K + 1;
                              e.Join( n1, n2 );    }    }    }

               // Get reps.

               vec<int> reps;
               e.OrbitRepsAlt(reps);
               int nr = reps.size( );

               // Build genomic locations.

               vec<vec<int>> rlocs(nr);
               for ( int i = 0; i < N; i++ )
               {    int p = BinPosition( reps, e.ClassId(i) );
                    for ( int j = 0; j < locs[i].isize( ); j++ )
                         rlocs[p].push_back( locs[i][j] );    }
               for ( int p = 0; p < nr; p++ )
                    UniqueSort( rlocs[p] );

               // Set up graph structure.

               vec<vec<int>> from(nr), to(nr);
               cout << Date( ) << ": building from and to" << endl;
               for ( int i = 0; i < ids.isize( ); i++ )
               {    int n = 0;
                    for ( int l = 0; l < i; l++ )
                         n += all[ ids[l] ].isize( ) - K + 1;
                    for ( int j = 0; j < all[ ids[i] ].isize( ) - K; j++ )
                    {    int x1 = BinPosition( reps, e.ClassId( n + j ) );
                         int x2 = BinPosition( reps, e.ClassId( n + j + 1 ) );
                         from[x1].push_back(x2), to[x2].push_back(x1);    }    }
               for ( int i = 0; i < nr; i++ )
               {    Sort( from[i] ), Sort( to[i] );    }

               // Make the graph.  Use as vertices pairs (kmer,locs) where locs
               // denote genomic locations.

               vec< pair<String,vec<int>> > greps( reps.size( ) );
               for ( int i = 0; i < reps.isize( ); i++ )
                    greps[i] = make_pair( kmers[ reps[i] ], rlocs[i] );
               digraphV< pair<String,vec<int>> > G( from, to, greps );

               // Zipper.

               while(1)
               {    Bool improved = False;
                    for ( int v = 0; v < G.N( ); v++ )
                    {    for ( int i1 = 0; i1 < G.From(v).isize( ); i1++ )
                         {    for ( int i2 = i1 + 1; i2 < G.From(v).isize( ); i2++ )
                              {    int w1 = G.From(v)[i1], w2 = G.From(v)[i2];
                                   if ( G.Vert(w1).first != G.Vert(w2).first ) 
                                        continue;
                                   if ( w1 == w2 ) continue;
                                   G.TransferEdges( w2, w1 );
                                   G.VertMutable(w1).second.append( 
                                        G.Vert(w2).second );
                                   UniqueSort( G.VertMutable(w1).second );
                                   cout << "zippered " << w1 << " to " << w2 
                                        << endl;    
                                   improved = True;
                                   goto next1;    }    }
                         next1:
                         for ( int i1 = 0; i1 < G.To(v).isize( ); i1++ )
                         {    for ( int i2 = i1 + 1; i2 < G.To(v).isize( ); i2++ )
                              {    int w1 = G.To(v)[i1], w2 = G.To(v)[i2];
                                   if ( G.Vert(w1).first != G.Vert(w2).first ) 
                                        continue;
                                   if ( w1 == w2 ) continue;
                                   G.TransferEdges( w2, w1 );
                                   G.VertMutable(w1).second.append( 
                                        G.Vert(w2).second );
                                   UniqueSort( G.VertMutable(w1).second );
                                   cout << "zippered " << w1 << " to " << w2 
                                        << endl;    
                                   improved = True;
                                   goto next2;    }    }
                         next2: continue;    }
                    if ( !improved ) break;
                    G.RemoveEdgelessVertices( );    }

               // Check for cycles.

               cout << Date( ) << ": checking for cycles" << endl;
               if ( G.Acyclic( ) ) cout << "acyclic" << endl;
               else cout << "has cycle" << endl;
               vec<int> core;
               G.CyclicCore(core);
               cout << "core = " << printSeq(core) << endl;
               cout << Date( ) << ": done" << endl;

               vec<String> gkmers;
               for ( int i = 0; i < G.N( ); i++ )
               {    String s = G.Vert(i).first;
                    for ( int j = 0; j < G.Vert(i).second.isize( ); j++ )
                         s += "\n" + ToString( G.Vert(i).second[j] );
                    gkmers.push_back(s);    }

               vec<int> vs;

               for ( int i = 0; i < G.N( ); i++ )
               {    if ( Member( G.Vert(i).second, 72849 ) )
                         vs.push_back(i);    }
               // vs.push_back( BinPosition( reps, e.ClassId(3000) ) );
               // vs.push_back( G.N( ) / 2 );

               int depth = 200;
               for ( int d = 0; d < depth; d++ )
               {    int nv = vs.size( );
                    for ( int i = 0; i < nv; i++ )
                    {    for ( int j = 0; j < G.From( vs[i] ).isize( ); j++ )
                              vs.push_back( G.From( vs[i] )[j] );
                         for ( int j = 0; j < G.To( vs[i] ).isize( ); j++ )
                              vs.push_back( G.To( vs[i] )[j] );    }
                    UniqueSort(vs);    }

               {    Ofstream( out, "/wga/dev/jaffe/BroadCRD/xxx.dot" );
                    // vec<String> color( G.N( ), "" );
                    // for ( int i = 0; i < core.isize( ); i++ )
                    //      color[ core[i] ] = "red";
                    // G.DOT_vl( out, gkmers, vs, color )

                    vec<String> color( G.N( ), "" );
                    for ( int i = 0; i < G.N( ); i++ )
                         if ( gkmers[i].Contains( "\n" ) ) color[i] = "red";
                    G.DOT_vl( out, gkmers, vs, color );

                    // G.DOT_vl( out, gkmers, vs );    
                         }    }
                    // G.DOT_vl( out, gkmers, "", vec<String>( ), "" );    }    }

          // Find repeated kmers.

          vec< triple<String,int,int> > pkmers;
          vec<int> ids;
          String s1 = all[id1].ToString( );
          for ( int p = 0; p <= all[id1].isize( ) - K; p++ )
          {    pkmers.push( s1.substr( p, K ), p, id1 );
               ids.push_back(id1);    }
          for ( int m = 0; m < alignsx[id1].isize( ); m++ )
          {    int id2 = alignsx[id1][m].first;
               String s2 = all[id2].ToString( );
               const align& a = alignsx[id1][m].second;
               int p1 = a.pos1( ), p2 = a.pos2( );
               vec< pair<int,int> > p2p1;
               for ( int j = 0; j < a.Nblocks( ); j++ ) 
               {    p2p1.push( p2, p1 );
                    if ( a.Gaps(j) > 0 )
                    {    for ( int k = 0; k < a.Gaps(j); k++ )
                         {    p2p1.push( p2, p1 );
                              p2++;    }    }
                    if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
                    p2p1.push( p2, p1 );
                    for ( int x = 0; x < a.Lengths(j); x++ )
                    {    ++p1;
                         ++p2;     
                         p2p1.push( p2, p1 );    }    }
               UniqueSort(p2p1);
               for ( int j = 0; j < p2p1.isize( ); j++ )
               {    int p2 = p2p1[j].first, p1 = p2p1[j].second;
                    int m;
                    for ( m = j + 1; m < p2p1.isize( ); m++ )
                         if ( p2p1[m].first != p2p1[j].first ) break;
                    if ( p2 > all[id2].isize( ) - K ) continue;
                    pkmers.push( s2.substr( p2, K ), p1, id2 );
                    ids.push_back(id2);
                    j = m - 1;    }    }
          UniqueSort(pkmers);
          vec<double> qpos;
          vec< triple< String, vec<int>, vec<int> > > qkmers;
          for ( int i = 0; i < pkmers.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < pkmers.isize( ); j++ )
               {    if ( pkmers[j].first != pkmers[i].first ) break;
                    if ( pkmers[j].second - pkmers[j-1].second > MAX_FUDGE ) 
                         break;    }
               vec<int> xids;
               for ( int k = i; k < j; k++ )
                    xids.push_back( pkmers[k].third );
               UniqueSort(xids);
               if ( xids.size( ) > 1 )
               {    vec<int> x;
                    for ( int k = i; k < j; k++ )
                         x.push_back( pkmers[k].second );
                    vec<int> y;
                    for ( int k = i; k < j; k++ )
                         y.push_back( pkmers[k].third );
                    qpos.push_back( Mean(x) );
                    qkmers.push( pkmers[i].first, x, y );    }
               i = j - 1;    }
          SortSync( qpos, qkmers );
          for ( int i = 0; i < qkmers.isize( ); i++ )
          {    basevector b( qkmers[i].first );
               kmer<K> x(b);

               vec<String> z;
               for ( int l = 0; l < qkmers[i].second.isize( ); l++ )
               {    z.push_back( ( qkmers[i].third[l] < nreads ? "+" : "-" )
                         + ToString( qkmers[i].second[l] ) );    }

               cout << qkmers[i].first << " " << printSeq(z);
               int64_t low = LowerBound1( gkmers_plus, x );
               int64_t high = UpperBound1( gkmers_plus, x );
               if ( loc_g[id1] >= 0 )
               {    int g = loc_g[id1];
                    int start = loc_start[id1], stop = loc_stop[id1];
                    cout << "   ";
                    for ( int64_t j = low; j < high; j++ )
                    {    if ( gkmers_plus[j].second == g
                              && gkmers_plus[j].third >= start
                              && gkmers_plus[j].third < stop )
                         {    cout << " " << gkmers_plus[j].second << "."
                                   << gkmers_plus[j].third;    }    }    }
               else cout << "   " << "!";
               cout << endl;    }    }

     cout << Date( ) << ": done" << endl;
     Scram(0);    }
