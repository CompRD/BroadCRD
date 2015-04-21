///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Try to capture the sequence within a gap.  Highly experimental, based on
// several test cases.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Intvector.h"
#include "MainTools.h"
#include "PrintAlignment.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/RemodelGapTools.h"
#include "paths/Unipath.h"
#include "paths/long/LongReadsToPaths.h"
#include "paths/long/ReadPath.h"
#include "paths/long/ReadStack.h"
#include "paths/long/large/GapToyTools.h"

// Delete weakly supported edges.

template<int K> void CleanGraph( const HyperBasevectorX& hb,
     const vecbasevector& bases, const vec<int>& mreach_ids,
     const int end_use, HyperBasevector& hbg, vec<int>& to_left, vec<int>& to_right )
{
     cout << Date( ) << ": deleting weakly supported edges" << endl;
     const int zpasses = 2;
     for ( int zpass = 1; zpass <= zpasses; zpass++ )
     {    vec<vec<int>> count0( hbg.E( ) );
          int left = -1;
          vec<int> rights( mreach_ids.size( ), -1 );
          vecbasevector woof;

          const basevector& E1 = hb.EdgeObject( mreach_ids[0] );
          int n1 = Min( E1.isize( ), end_use );
          basevector b1( E1, E1.isize( ) - n1, n1 );
          woof.push_back(b1);
          for ( int i = 1; i < mreach_ids.isize( ); i++ )
          {    const basevector& E = hb.EdgeObject( mreach_ids[i] );
               basevector b( E, 0, Min( E.isize( ), end_use ) );
               woof.push_back(b);    }

          woof.Append(bases);
          for ( int e = 0; e < hbg.E( ); e++ )
               woof.push_back( hbg.EdgeObject(e) );
          vec< triple<kmer<K>,int,int> > kmers_plus;
          MakeKmerLookup0( woof, kmers_plus );
          for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
          {    int64_t j;
               for ( j = i + 1; j < kmers_plus.jsize( ); j++ )
                    if ( kmers_plus[j].first != kmers_plus[i].first ) break;
               for ( int64_t k = i; k < j; k++ )
               {    int e = kmers_plus[k].second, epos = kmers_plus[k].third;
                    if ( e == 0 && epos == 0 )
                    {    for ( int64_t l = i; l < j; l++ )
                         {    int e2 = kmers_plus[l].second 
                                   - mreach_ids.isize( ) - (int) bases.size( );
                              if ( e2 >= 0 ) left = to_left[e2];    }    }
                    if ( e > 0 && e < mreach_ids.isize( ) 
                         && epos == woof[e].isize( ) - K )
                    {    for ( int64_t l = i; l < j; l++ )
                         {    int e2 = kmers_plus[l].second 
                                   - mreach_ids.isize( ) - (int) bases.size( );
                              if ( e2 >= 0 ) rights[e] = to_right[e2];    }    }    }
               for ( int64_t k1 = i; k1 < j; k1++ )
               for ( int64_t k2 = i; k2 < j; k2++ )
               {    int id1 = kmers_plus[k1].second - mreach_ids.isize( ); // read id
                    if ( id1 < 0 ) continue;
                    if ( id1 >= (int) bases.size( ) ) continue;
                    int id2 = kmers_plus[k2].second 
                         - mreach_ids.isize( ) - (int) bases.size( );
                    if ( id2 < 0 ) continue;
                    count0[id2].push_back(id1);    }
               i = j - 1;    }
          for ( int e = 0; e < hbg.E( ); e++ )
               UniqueSort( count0[e] );
          vec<int> dels;
          if ( zpass >= 2 )
          {
          for ( int v = 0; v < hbg.N( ); v++ )
          {    if ( hbg.From(v).size( ) != 2 ) continue;
               int e1 = hbg.IFrom( v, 0 ), e2 = hbg.IFrom( v, 1 );
               if ( count0[e1].isize( ) < count0[e2].isize( ) ) swap( e1, e2 );
               if ( count0[e1].isize( ) >= 10 && count0[e2].isize( ) <= 1 ) 
                    dels.push_back(e2);
               if ( count0[e1].isize( ) >= 5 && count0[e2].isize( ) <= 0 ) 
                    dels.push_back(e2);    }
          for ( int v = 0; v < hbg.N( ); v++ )
          {    if ( hbg.To(v).size( ) != 2 ) continue;
               int e1 = hbg.ITo( v, 0 ), e2 = hbg.ITo( v, 1 );
               if ( count0[e1].isize( ) < count0[e2].isize( ) ) swap( e1, e2 );
               if ( count0[e1].isize( ) >= 10 && count0[e2].isize( ) <= 1 ) 
                    dels.push_back(e2);
               if ( count0[e1].isize( ) >= 5 && count0[e2].isize( ) <= 0 ) 
                    dels.push_back(e2);    }
          }
          PRINT2( left, rights[1] );

          if ( left >= 0 )
          {    vec<int> bet;
               Bool have_right = False;
               for ( int ri = 1; ri < rights.isize( ); ri++ )
               {    if ( rights[ri] >= 0 )
                    {    have_right = True;
                         bet.append( hbg.EdgesSomewhereBetween( left, rights[ri] ) );
                              }    }
               UniqueSort(bet);
               if (have_right)
               for ( int e = 0; e < hbg.E( ); e++ )
                    if ( !BinMember( bet, e ) ) dels.push_back(e);    }

          hbg.DeleteEdges(dels);
          hbg.RemoveUnneededVertices( );
          hbg.RemoveDeadEdgeObjects( );
          PRINT( hbg.E( ) );
          hbg.ToLeft(to_left), hbg.ToRight(to_right);    }    }

int main( int argc, char *argv[] ) {
    RunTime( );

    BeginCommandArguments;
    CommandArgument_Int(e1);
    CommandArgument_Int(e2);
    EndCommandArguments;

     // Say what we're doing.

     cout << "Looking for path from e1 to e2 in 51400.newchem/a.final." << endl;
     cout << "Will generate yyy.{dot,fasta} in ~/crd." << endl << endl;

     // Control.

     Bool print_walks = False; // reads versus graph paths

     // Define directory.

     // String work_dir = "/wga/scr4/jaffe/GapToy/51400.newchem";
     // String work_dir = "/wga/scr4/jaffe/GapToy/52009.HCC1143+BL";
     String work_dir = "/wga/scr4/jaffe/GapToy/52398.XDP2";

     // Constants.

     const int read_len = 250; // *************************************************
     const int end_back = 400;
     const int K = 100;
     const int M = 28;
     const int L = 40;
     const int end_use = 1000;

     vec<int> passes = {0,1}; // {a.200,a.final} -- change when examples added

     int mpasses = passes.size( );
     int total_passes = Sum(passes);
     int pass_id = 0;
     const int mpass = 2;
     {    String as = "a.final";
          String dir = work_dir + "/" + as;
          // ---->

     // Load assembly.

     cout << Date( ) << ": loading assembly" << endl;
     HyperBasevectorX hb;
     BinaryReader::readFile( dir + "/a.hbx", &hb );
     vec<int> inv;
     BinaryReader::readFile( dir + "/a.inv", &inv );
     
     // Passes for examples.

     for ( int pass = 1; pass <= passes[mpass-1]; pass++ )
     {

     ++pass_id;
     cout << Date( ) << ": start pass " << pass << endl;
     PRINT(pass_id);

     // Pick an a.final edge that dead ends on the right.

     basevector ref;
     int nedges = -1;

     int did = e1;

     // Identify the tenx barcodes, then load the reads.

     vec<uint32_t> bx;
     vec<vec<int64_t>> hits;
     vec<int64_t> all;
     vecbasevector tenx_bases;

     // Find the ids of the reads that are incident upon the dead end.

     cout << Date( ) << ": finding read ids" << endl;
     vec<int> s = { did, inv[did] };
     VecULongVec x;
     x.Read( dir + "/a.paths.inv", s );
     vec<int64_t> xall;
     for ( int i = 0; i < (int) x.size( ); i++ )
     for ( int j = 0; j < (int) x[i].size( ); j++ )
     {    int64_t pid = x[i][j]/2;
          xall.push_back( 2*pid, 2*pid+1 );    }
     UniqueSort(xall);

     // Load the paths for these reads.

     cout << Date( ) << ": loading " << xall.size( ) << " paths" << endl;
     ReadPathVec paths, paths2;
     paths.Read( dir + "/a.paths", xall );

     // Determine what the reads that reach off the end reach to.

     cout << Date( ) << ": determine where we reach" << endl;
     vec< pair<int,int> > reach; // (edge, start rel. end of edge did)
     int frag_size = 455; // rough, could compute
     for ( int i = 0; i < (int) paths.size( ); i++ )
     {    const ReadPath& p = paths[i];
          if ( p.size( ) == 0 ) continue;
          int ip = ( i % 2 == 0 ? i+1 : i-1 );
          int start = p.getOffset( );
          for ( int j = 0; j < (int) p.size( ); j++ )
          {    int e = p[j];
               if ( e == did )
               {    int stop = start + read_len;
                    if ( stop > hb.Bases(e) - end_back )
                    {    if ( paths[ip].size( ) > 0 )
                         {    int startp = paths[ip].getOffset( );
                              for ( int l = 0; l < (int) paths[ip].size( ); l++ )
                              {    int ep = inv[ paths[ip][l] ];
                                   if ( ep == e ) continue;
                                   int gap = frag_size - ( hb.Bases(e) - start )
                                        - ( hb.Bases(ep) - startp );
                                   reach.push( ep, gap );
                                   startp -= hb.Kmers(ep);    }    }    }    }
               start -= hb.Kmers(e);    }    }
     Sort(reach);
     vec< pair<int,int> > mreach;
     for ( int i = 0; i < reach.isize( ); i++ )
     {    int j, sum = 0;
          for ( j = i + 1; j < reach.isize( ); j++ )
               if ( reach[j].first != reach[i].first ) break;
          for ( int k = i; k < j; k++ )
               sum += reach[k].second;
          int start = int( round( double(sum)/(j-i) ) );
          cout << "a) " << j-i << " links to " << reach[i].first << endl; // XXXXXXX
          mreach.push( reach[i].first, start );
          i = j - 1;    }
     mreach.push( e2, 0 );
     Bool have_did = False;
     for ( int i = 0; i < mreach.isize( ); i++ )
          if ( mreach[i].first == did ) have_did = True;
     if ( !have_did ) mreach.push( did, -hb.Bases(did) );
     int mp = -1;
     for ( int i = 0; i < mreach.isize( ); i++ )
          if ( mreach[i].first == did ) mp = i;
     if ( mp > 0 ) swap( mreach[mp], mreach[0] );

     // Summarize.

     cout << "\nstarting from edge " << did << endl;
     cout << "edge, start for edges in a.200\n";
     for ( int i = 0; i < mreach.isize( ); i++ )
          cout << mreach[i].first << ", " << mreach[i].second << endl;

     // Find the ids of all the reads on any of the edges.

     cout << Date( ) << ": finding all the read ids" << endl;
     vec<int64_t> xalle;
     {    vec<int> s;
          for ( int i = 0; i < mreach.isize( ); i++ )
               s.push_back( mreach[i].first, inv[ mreach[i].first ] );
          UniqueSort(s);
          VecULongVec x;
          x.Read( dir + "/a.paths.inv", s );
          vec<int64_t> xall;
          for ( int i = 0; i < (int) x.size( ); i++ )
          for ( int j = 0; j < (int) x[i].size( ); j++ )
          {    int64_t pid = x[i][j]/2;
               xalle.push_back( 2*pid, 2*pid+1 );    }    }
     UniqueSort(xalle);

     // Now load the paths for these.

     cout << Date( ) << ": next loading " << xalle.size( ) << " paths" << endl;
     ReadPathVec pathse;
     pathse.Read( dir + "/a.paths", xalle );
     
     // Now gather up all the reads, and compute (id, orientation).

     cout << Date( ) << ": looking for gappers" << endl;
     vec< pair<int64_t,Bool> > gappers;
     for ( int i = 0; i < xalle.isize( ); i++ )
     {    const ReadPath& p = pathse[i];
          if ( p.size( ) == 0 ) continue;
          int64_t id1 = xalle[i];
          int64_t id2 = ( id1 % 2 == 0 ? id1 + 1 : id1 - 1 );
          int start = p.getOffset( );
          for ( int j = 0; j < (int) p.size( ); j++ )
          {    int e = p[j];
               for ( int l = 0; l < mreach.isize( ); l++ )
               {    if ( mreach[l].first == e )
                    {    int mstart = mreach[l].second;
                         int rstart = start + mstart;
                         if ( rstart + 250 > - end_back && rstart < end_back )
                         {    gappers.push( id1, True );
                              gappers.push( id2, False );    }    }
                    if ( mreach[l].first == inv[e] )
                    {    int mstart = mreach[l].second;
                         int rstart = mstart + ( hb.Bases(e) - start - 250 );
                         if ( rstart > -400 && rstart < 700 )
                         {    gappers.push( id1, False );
                              gappers.push( id2, True );    }    }    }
               start -= hb.Kmers(e);    }    }
     UniqueSort(gappers);

     // Print the gappers.

     Bool print_gappers = False;
     if (print_gappers)
     {    for ( int i = 0; i < gappers.isize( ); i++ )
          {    int64_t id = gappers[i].first;
               cout << ( gappers[i].second ? "+" : "-" ) << gappers[i].first;
               int p = BinPosition( xalle, id );
               if ( p >= 0 && pathse[p].size( ) > 0 )
               {    cout << " a.final=" << pathse[p].getOffset( )
                         << ":" << printSeq(pathse[p]);    }
               cout << endl;    }    }

     // Load and orient the gappers.

     cout << Date( ) << ": loading " << gappers.size( ) << " reads" << endl;
     vecbasevector bases;
     VecPQVec qualsp;
     vec<int64_t> gids( gappers.size( ) );
     for ( int i = 0; i < (int) gappers.size( ); i++ )
          gids[i] = gappers[i].first;
     bases.Read( dir + "/../data/frag_reads_orig.fastb", gids );
     qualsp.Read( dir + "/../data/frag_reads_orig.qualp", gids );
     vecqualvector quals( qualsp.size( ) );
     for ( int i = 0; i < (int) quals.size( ); i++ )
          qualsp[i].unpack( &quals[i] );
     for ( int i = 0; i < (int) gappers.size( ); i++ )
     {    if ( !gappers[i].second ) 
          {    bases[i].ReverseComplement( );
               quals[i].ReverseMe( );    }    }
     bases.WriteAll( "fido" + ToString(pass_id) + ".fastb" );
     quals.WriteAll( "fido" + ToString(pass_id) + ".qualb" );

     // Align the gappers to each other.

     cout << Date( ) << ": aligning" << endl;
     vec< triple<int,int,int> > aligns;
     vec< triple<kmer<L>,int,int> > kmers_plus;
     MakeKmerLookup0( bases, kmers_plus );
     for ( int i = 0; i < kmers_plus.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < kmers_plus.isize( ); j++ )
               if ( kmers_plus[j].first != kmers_plus[i].first ) break;
          for ( int k1 = i; k1 < j; k1++ )
          for ( int k2 = i; k2 < j; k2++ )
          {    if ( k2 == k1 ) continue;
               int id1 = kmers_plus[k1].second, id2 = kmers_plus[k2].second;
               int offset = kmers_plus[k1].third - kmers_plus[k2].third;
               aligns.push( id1, id2, offset );    }
          i = j - 1;    }
     UniqueSort(aligns);
     PRINT( aligns.size( ) );

     // Remove aligns having hq mismatches.

     vec<Bool> to_delete( aligns.size( ), False );
     for ( int i = 0; i < aligns.isize( ); i++ )
     {    vec<int> qdiffs;
          int id1 = aligns[i].first, id2 = aligns[i].second;
          int offset = aligns[i].third;
          const basevector &b1 = bases[id1], &b2 = bases[id2];
          const qualvector &q1 = quals[id1], &q2 = quals[id2];
          if ( offset >= 0 )
          {    for ( int j = 0; j < b2.isize( ); j++ )
               {    if ( j+offset >= b1.isize( ) ) break;
                    if ( b2[j] != b1[j+offset] )
                         qdiffs.push_back( Min( q2[j], q1[j+offset] ) );    }    }
          else
          {    for ( int j = 0; j < b1.isize( ); j++ )
               {    if ( j-offset >= b2.isize( ) ) break;
                    if ( b1[j] != b2[j-offset] )
                         qdiffs.push_back( Min( q1[j], q2[j-offset] ) );    }    }
          ReverseSort(qdiffs);
          int qmax = ( qdiffs.empty( ) ? 0 : Max(qdiffs) );
          if ( qmax >= 30 ) to_delete[i] = True;    }
     EraseIf( aligns, to_delete );
     PRINT(aligns.size( ) );

     // Print alignments.

     Bool print_aligns = False;
     if (print_aligns)
     {    for ( int i = 0; i < aligns.isize( ); i++ )
          {    int id1 = aligns[i].first, id2 = aligns[i].second;
               cout << "\nalignment of reads " << id1 << " and " << id2 << endl;
               int offset = aligns[i].third;
               const basevector &b1 = bases[id1], &b2 = bases[id2];
               const qualvector &q1 = quals[id1], &q2 = quals[id2];
               int pos1, pos2;
               if ( offset >= 0 )
               {    pos1 = offset;
                    pos2 = 0;    }
               else
               {    pos1 = 0;
                    pos2 = -offset;    }
               avector<int> gaps(1), lengths(1);
               gaps(0) = 0;
               lengths(0) = IntervalOverlap( 
                    0, b1.isize( ), offset, offset + b2.isize( ) );
               align a( pos1, pos2, gaps, lengths );
               PrintVisualAlignment( True, cout, b1, b2, a, q1, q2 );    }    }

     // Build stack for each read, and correct it.

     vecbasevector basesx( bases.size( ) );
     vecqualvector qualsx( quals.size( ) );
     vec<readstack> stacks( bases.size( ) );
     for ( int id1 = 0; id1 < (int) bases.size( ); id1++ )
     {    
          // Let x = { ( id, offset ) }.

          vec< pair<int,int> > x;
          for ( int j = 0; j < aligns.isize( ); j++ )
          {    if ( aligns[j].first != id1 ) continue;
               x.push( aligns[j].third, aligns[j].second );    }
          Sort(x);
          x.push_front( make_pair( 0, id1 ) );

          // If the read is +, build stack extending to the right.
          // If the read is -, build stack extending to the left.

          readstack stack;
          int n = x.size( ), k = bases[id1].size( );
          if ( gappers[id1].second )
          {    for ( int j = 0; j < x.isize( ); j++ )
                    k = Max( k, x[j].first + bases[ x[j].second ].isize( ) );
               stack.Initialize( n, k );
               for ( int j = 0; j < n; j++ )
               {    int id2 = x[j].second, offset = x[j].first;
                    const basevector& b2 = bases[id2];
                    const qualvector& q2 = quals[id2];
                    for ( int p2 = 0; p2 < b2.isize( ); p2++ )
                    {    int p1 = p2 + offset;
                         if ( p1 < 0 ) continue;
                         stack.SetBase( j, p1, b2[p2] );
                         stack.SetQual( j, p1, q2[p2] );     
                         stack.SetOffset( j, offset );
                         stack.SetLen( j, b2.size( ) );
                         stack.SetId( j, id2 );
                         stack.SetRc2( j, False );    }    }    }
          else
          {    int left = 0;
               for ( int j = 0; j < x.isize( ); j++ )
                    left = Min( left, x[j].first );
               int k = bases[id1].isize( ) - left;
               stack.Initialize( n, k );
               for ( int j = 0; j < n; j++ )
               {    int id2 = x[j].second, offset = x[j].first;
                    const basevector& b2 = bases[id2];
                    const qualvector& q2 = quals[id2];
                    for ( int p2 = 0; p2 < b2.isize( ); p2++ )
                    {    int p1 = p2 + offset;
                         if ( p1 - left >= k ) continue;
                         stack.SetBase( j, p1 - left, b2[p2] );
                         stack.SetQual( j, p1 - left, q2[p2] );     
                         stack.SetOffset( j, offset );
                         stack.SetLen( j, b2.size( ) );
                         stack.SetId( j, id2 );
                         stack.SetRc2( j, False );    }    }    }
          stacks[id1] = stack;

          // Get consensus.
          
          // if ( !gappers[id1].second ) stack.Reverse( );
          Bool print_stack = False;
          if (print_stack)
          {    cout << "\nstack for read " << id1 << endl;
               stack.Print(cout);    }
          int trim_to;
          // stack.CorrectAll( basesx[id1], qualsx[id1], trim_to );
          stack.Consensus1( basesx[id1], qualsx[id1] );
               }

     // Try to merge each pair.

     vec< vec< pair<int,int> > > offsets( bases.size( ) / 2 );
     vecbasevector frags;
     Bool print_basic = False;
     for ( int pid = 0; pid < (int) bases.size( ) / 2; pid++ )
     {    if (print_basic)
               cout << "\npair " << pid << " = " << gappers[2*pid].first/2 << endl;
          int id1 = 2*pid, id2 = 2*pid+1;
          if ( !gappers[id1].second ) swap( id1, id2 );
          int offset = -1;
          vec<int> off = GetOffsets1( stacks[id1], stacks[id2], 0, 0 );
          if ( off.solo( ) )
          {    offset = off[0];
               if (print_basic)
               {    cout << "using offset " << offset << endl;
                    cout << "predict offsets = " << printSeq(off) << endl;
                    cout << "fragment size = " << offset + basesx[id2].isize( ) 
                         << endl;
                    cout << ( gappers[id1].second ? "fw" : "rc" ) << " : "
                         << ( gappers[id2].second ? "fw" : "rc" ) << endl;    }    }
     
          Bool print_reads = False;
          if (print_reads)
          {    bases[id1].Print( 
                    cout, ToString(pid) + ".1.orig" + " = " + ToString(id1) );
               bases[id2].Print( 
                    cout, ToString(pid) + ".2.orig" + " = " + ToString(id2) );
               basesx[id1].Print( 
                    cout, ToString(pid) + ".1" + " = " + ToString(id1) );
               basesx[id2].Print( 
                    cout, ToString(pid) + ".2" + " = " + ToString(id2) );    }
          if ( !off.solo( ) ) continue;

          readstack s( stacks[id1] );
          s.Merge( stacks[id2], offset );
          basevector f;
          qualvector fragq;
          s.Consensus1( f, fragq );

          if (print_reads)
          {    cout << "found fragment size = " << f.size( ) << endl;
               f.Print( cout, ToString(pid) + ".join" );    }
          frags.push_back(f);    }

     // Form the graph from the two edges and the reads.

     cout << Date( ) << ": forming graph" << endl;
     HyperBasevector hbg;

     {    vecbasevector stuff;

          const basevector& E1 = hb.EdgeObject( mreach[0].first );
          int n1 = Min( E1.isize( ), end_use );
          basevector b1( E1, E1.isize( ) - n1, n1 );
          stuff.push_back(b1);
          for ( int i = 1; i < mreach.isize( ); i++ )
          {    const basevector& E = hb.EdgeObject( mreach[i].first );
               basevector b( E, 0, Min( E.isize( ), end_use ) );
               stuff.push_back(b);    }

          // for ( int i = 0; i < mreach.isize( ); i++ )
          //      stuff.push_back( hb.EdgeObject( mreach[i].first ) );

          stuff.Append(frags);
          vecKmerPath paths, paths_rc, unipaths;
          // vec<big_tagged_rpint> pathsdb, unipathsdb;
          vec<tagged_rpint> pathsdb, unipathsdb;
          ReadsToPathsCoreY( stuff, K, paths );
          CreateDatabase( paths, paths_rc, pathsdb );
          Unipath( paths, paths_rc, pathsdb, unipaths, unipathsdb );
          digraph A;
          BuildUnipathAdjacencyGraph( paths, paths_rc, pathsdb, unipaths,
               unipathsdb, A );
          HyperKmerPath h;
          BuildUnipathAdjacencyHyperKmerPath( K, A, unipaths, h );
          KmerBaseBroker kbb( K, paths, paths_rc, pathsdb, stuff );
          hbg = HyperBasevector( h, kbb );    }
     vec<int> to_left, to_right;
     hbg.ToLeft(to_left), hbg.ToRight(to_right);

     // Delete weakly supported edges.

     cout << Date( ) << ": deleting weakly supported edges" << endl;
     vec<int> mreach_ids;
     for ( int i = 0; i < mreach.isize( ); i++ )
          mreach_ids.push_back( mreach[i].first );
     CleanGraph<K>( hb, bases, mreach_ids, end_use, hbg, to_left, to_right );

     // Now find the paths that each read pair could represent in the graph.

     vec<int> gat( gappers.size( ), -1 );
     vec<int> gatpos( gappers.size( ), -1 );
     {    vecbasevector woof;
          woof.Append(bases);
          for ( int e = 0; e < hbg.E( ); e++ )
               woof.push_back( hbg.EdgeObject(e) );
          vec< triple<kmer<K>,int,int> > kmers_plus;
          MakeKmerLookup0( woof, kmers_plus );
          for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
          {    int64_t j;
               for ( j = i + 1; j < kmers_plus.jsize( ); j++ )
                    if ( kmers_plus[j].first != kmers_plus[i].first ) break;
               for ( int64_t k1 = i; k1 < j; k1++ )
               for ( int64_t k2 = i; k2 < j; k2++ )
               {    int id1 = kmers_plus[k1].second; // read id
                    if ( id1 < 0 ) continue;
                    if ( id1 >= (int) bases.size( ) ) continue;
                    int pos1 = kmers_plus[k1].third;
                    Bool fw1 = gappers[id1].second;
                    if ( fw1 && pos1 != 0 ) continue;
                    if ( !fw1 && pos1 != bases[id1].isize( ) - K ) continue;
                    int id2 = kmers_plus[k2].second - (int) bases.size( );
                    if ( id2 < 0 ) continue;
                    int pos2 = kmers_plus[k2].third;
                    gat[id1] = id2;
                    if (fw1) gatpos[id1] = pos2;
                    else gatpos[id1] = pos2 + K;
                    // count0[id2].push_back(id1);    
                         }
               i = j - 1;    }    }

     // Walk the reads.

     cout << "\nwalking reads" << endl;
     HyperBasevector hbgr(hbg);
     hbgr.Reverse( );
     vec<vec<int>> all_paths( gappers.size( ) );
     for ( int id = 0; id < gappers.isize( ); id++ )
     {    if ( gat[id] < 0 ) continue;
          vec<int> e = { gat[id] };
          const int min_gain = 5;
          const int mode = 1;
          int offset;
          basevector b;
          vec<int> pv;
          basevector b1 = bases[id];
          qualvector q1 = quals[id];
          ReadPath p;

          if ( gappers[id].second )
          {    p = ReadPath( gatpos[id], e );
               ExtendPath( p, id, hbg, to_right, b1, q1, min_gain, False, mode );
               if (print_walks)
               {    cout << "+" << id << " --> " << p.getOffset( ) << ": " 
                         << printSeq(p) << endl;     }
               for ( int i = 0; i < (int) p.size( ); i++ )
                    pv.push_back( p[i] );
               b = hbg.Cat(pv);
               offset = -p.getOffset( );    }

          else
          {    p = ReadPath( hbg.Bases( e[0] ) - gatpos[id], e );
               b1.ReverseComplement( );
               q1.ReverseMe( );
               ExtendPath( p, id, hbgr, to_left, b1, q1, min_gain, False, mode );
               if (print_walks)
               {    cout << "-" << id << " --> " << p.getOffset( ) << ": " 
                         << printSeq(p) << endl;     }
               for ( int i = 0; i < (int) p.size( ); i++ )
                    pv.push_back( p[i] );
               b = hbgr.Cat(pv);
               offset = -p.getOffset( );    }

          int pos1, pos2;
          if ( offset >= 0 )
          {    pos1 = offset;
               pos2 = 0;    }
          else
          {    pos1 = 0;
               pos2 = -offset;    }
          const basevector &b2 = b;
          avector<int> gaps(1), lengths(1);
          gaps(0) = 0;
          lengths(0) = IntervalOverlap(
               0, b1.isize( ), offset, offset + b2.isize( ) );
          align a( pos1, pos2, gaps, lengths );
          vec<int> qdiffs;
          for ( int l = 0; l < lengths(0); l++ )
               if ( b1[pos1+l] != b[pos2+l] ) qdiffs.push_back( q1[pos1+l] );
          ReverseSort(qdiffs);
          if (print_walks) cout << "qdiffs: " << printSeq(qdiffs) << endl;
          if ( qdiffs.empty( ) || qdiffs[0] < 30 )
          {    vec<int> x;
               for ( int j = 0; j < (int) p.size( ); j++ )
                    x.push_back( p[j] );
               if ( !gappers[id].second ) x.ReverseMe( );
               all_paths[id] = x;    }
          if (print_walks) PrintVisualAlignment( True, cout, b1, b, a, q1 );    }

     // Summarize paths.

     vec<vec<int>> allp, allpu;
     for ( int i = 0; i < gappers.isize( ); i += 2 )
     {    if ( all_paths[i].nonempty( ) && all_paths[i+1].nonempty( ) )
               allp.push_back( all_paths[i], all_paths[i+1] );    }
     Sort(allp);
     cout << "\npath summary:\n";
     for ( int i = 0; i < allp.isize( ); i++ )
     {    int j = allp.NextDiff(i);
          allpu.push_back( allp[i] );
          cout << "[" << j-i << "] " << printSeq( allp[i] ) << endl;
          i = j - 1;    }
     cout << "\n";

     // Rebuild using paths.

     cout << Date( ) << ": rebuilding graph using large K" << endl;
     vec<vec<int>> fp;
     for ( int i = 0; i < allpu.isize( ); i++ )
     {    vec<int> p = allpu[i];
          int v = to_left[ allpu[i].front( ) ], w = to_right[ allpu[i].back( ) ];
          if ( hbg.To(v).solo( ) ) p.push_front( hbg.ITo( v, 0 ) );
          if ( hbg.From(w).solo( ) ) p.push_back( hbg.IFrom( w, 0 ) );
          fp.push_back(p);    }
     UniqueSort(fp);
     {    vecbasevector stuff;

          const basevector& E1 = hb.EdgeObject( mreach[0].first );
          int n1 = Min( E1.isize( ), end_use );
          basevector b1( E1, E1.isize( ) - n1, n1 );
          stuff.push_back(b1);
          for ( int i = 1; i < mreach.isize( ); i++ )
          {    const basevector& E = hb.EdgeObject( mreach[i].first );
               basevector b( E, 0, Min( E.isize( ), end_use ) );
               stuff.push_back(b);    }

          for ( int i = 0; i < fp.isize( ); i++ )
               stuff.push_back( hbg.Cat( fp[i] ) );
          vecKmerPath paths, paths_rc, unipaths;
          vec<tagged_rpint> pathsdb, unipathsdb;
          const int K = 200;
          ReadsToPathsCoreY( stuff, K, paths );
          CreateDatabase( paths, paths_rc, pathsdb );
          Unipath( paths, paths_rc, pathsdb, unipaths, unipathsdb );
          digraph A;
          BuildUnipathAdjacencyGraph( paths, paths_rc, pathsdb, unipaths,
               unipathsdb, A );
          HyperKmerPath h;
          BuildUnipathAdjacencyHyperKmerPath( K, A, unipaths, h );
          KmerBaseBroker kbb( K, paths, paths_rc, pathsdb, stuff );
          hbg = HyperBasevector( h, kbb );
          hbg.ToLeft(to_left);
          hbg.ToRight(to_right);    }

     // Delete weakly supported edges. (AGAIN!)

     cout << Date( ) << ": deleting weakly supported edges" << endl;
     CleanGraph<200>( hb, bases, mreach_ids, end_use, hbg, to_left, to_right );

     // Test for reference.

     if ( ref.size( ) > 0 )
     {    cout << "\nvalidating" << endl;
          vec<int> sources, sinks;
          hbg.Sources(sources);
          hbg.Sinks(sinks);
          Bool have_ref = False;
          for ( int i1 = 0; i1 < sources.isize( ); i1++ )
          for ( int i2 = 0; i2 < sinks.isize( ); i2++ )
          {    vec<vec<int>> paths;
               hbg.EdgePaths( sources[i1], sinks[i2], paths );
               for ( int i = 0; i < paths.isize( ); i++ )
               {    basevector b = hbg.Cat( paths[i] );
                    // if ( b == ref )
                    if ( search( b.begin( ), b.end( ), ref.begin( ), ref.end( ) )
                         != b.end( ) )
                    {    cout << "Cool, found reference!" << endl;
                         have_ref = True;    }
                    if ( have_ref && hbg.E( ) <= nedges )
                    {    cout << pass_id << " VALIDATED!" << endl;    
                         goto validation_complete;    }    }    }    }
     validation_complete:

     // Identify edges sharing kmers with e1 and e2.

     cout << Date( ) << ": mapping back" << endl;
     vec<Bool> m1( hbg.E( ), False ), m2( hbg.E( ), False );
     vec<vec<int>> count( hbg.E( ) );
     {    vecbasevector woof;
          woof.push_back( hb.EdgeObject(e1) );
          woof.push_back( hb.EdgeObject(e2) );
          woof.Append(bases);
          for ( int e = 0; e < hbg.E( ); e++ )
               woof.push_back( hbg.EdgeObject(e) );
          vec< triple<kmer<K>,int,int> > kmers_plus;
          MakeKmerLookup0( woof, kmers_plus );
          for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
          {    int64_t j;
               for ( j = i + 1; j < kmers_plus.jsize( ); j++ )
                    if ( kmers_plus[j].first != kmers_plus[i].first ) break;

               for ( int64_t k1 = i; k1 < j; k1++ )
               for ( int64_t k2 = i; k2 < j; k2++ )
               {    int id1 = kmers_plus[k1].second;
                    if ( id1 >= mreach_ids.isize( ) ) continue;
                    int id2 = kmers_plus[k2].second - 2 - (int) bases.size( );
                    if ( id2 < 0 ) continue;
                    if ( id1 == 0 ) m1[id2] = True;
                    else m2[id2] = True;    }

               for ( int64_t k1 = i; k1 < j; k1++ )
               for ( int64_t k2 = i; k2 < j; k2++ )
               {    int id1 = kmers_plus[k1].second - 2; // read id
                    if ( id1 < 0 ) continue;
                    if ( id1 >= (int) bases.size( ) ) continue;
                    int id2 = kmers_plus[k2].second 
                         - 2 - (int) bases.size( );
                    if ( id2 < 0 ) continue;
                    count[id2].push_back(id1);    }

               i = j - 1;    }    }
     for ( int e = 0; e < hbg.E( ); e++ )
          UniqueSort( count[e] );

     // Dump fasta file.

     {    Ofstream( fout, "/wga/dev/jaffe/BroadCRD/yyy.fasta" );
          for ( int e = 0; e < hbg.E( ); e++ )
               hbg.EdgeObject(e).Print( fout, e );    }

     // Display the graph.

     cout << Date( ) << ": displaying graph" << endl;
     vec<String> edge_color( hbg.E( ) );
     for ( int e = 0; e < hbg.E( ); e++ )
     {    if ( m1[e] && m2[e] ) edge_color[e] = "purple";
          else if ( m1[e] ) edge_color[e] = "red";
          else if ( m2[e] ) edge_color[e] = "blue";    
               }
     vec<String> edge_names( hbg.E( ) );
     for ( int e = 0; e < hbg.E( ); e++ )
          edge_names[e] = ToString(e) + "(" + ToString( count[e].size( ) ) + ")";
     vec<double> lengths( hbg.E( ) );
     for ( int i = 0; i < hbg.E( ); i++ )
          lengths[i] = hbg.Kmers(i);
     {    Ofstream( dout, "/wga/dev/jaffe/BroadCRD/yyy.dot" );
          const Bool LABEL_CONTIGS = False;
          const Bool DOT_LABEL_VERTICES = False;
          hbg.PrettyDOT( dout, lengths, HyperBasevector::edge_label_info(
               HyperBasevector::edge_label_info::DIRECT, &edge_names ),
               LABEL_CONTIGS, DOT_LABEL_VERTICES, NULL, NULL, NULL, NULL,
               NULL, &edge_color );    }
     }    }

     cout << Date( ) << ": done" << endl;
     Scram(0);    }
