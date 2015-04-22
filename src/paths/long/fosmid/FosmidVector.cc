///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "CoreTools.h"
#include "FastIfstream.h"
#include "FetchReads.h"
#include "PairsManager.h"
#include "Qualvector.h"
#include "TokenizeString.h"
#include "VecUtilities.h"
#include "efasta/EfastaTools.h"
#include "lookup/LookAlign.h"
#include "paths/long/fosmid/FosmidPool.h"
#include "paths/long/fosmid/FosmidVector.h"
#include "paths/long/Logging.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/SupportedHyperBasevector.h"

namespace { // open anonymous namespace

// What follows are problem the Fosmid vector boundaries for hpool1.

basevector vleft(  "TCAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCA" );
basevector vright( "ACTATAGGGCGAATTCGAGCTCGGTACCCGGGGATCCCAC" );

} // end of anonymous namespace

void RemoveFinished( const String& SAMPLE, const String& TMP )
{    String fdir = "/wga/dev/references/Homo_sapiens/NA12878_Fosmid_Pool.regions.fin";
     vec<String> all = AllFiles(fdir);
     vecbasevector genome;
     for ( int i = 0; i < all.isize( ); i++ )
     {    if ( !all[i].Contains( ".fasta", -1 ) ) continue;
          int id = all[i].Between( ".", "." ).Int( );
          if ( SAMPLE == "hpool2" && id > 55 ) continue;
          if ( SAMPLE == "hpool3" && id <= 55 ) continue;
          vecbasevector b;
          FetchReads( b, 0, fdir + "/" + all[i] );
          genome.Append(b);
          for ( int j = 0; j < (int) b.size( ); j++ )
               genome.push_back( ReverseComplement( b[j] ) );    }

     vecbasevector bases( TMP + "/frag_reads_orig.fastb" );
     vecqualvector quals( TMP + "/frag_reads_orig.qualb" );
     PairsManager pairs;
     pairs.Read( TMP + "/frag_reads_orig.pairs" );

     const int K = 80;
     vec< triple<kmer<K>,int,int> > kmers_plus;
     MakeKmerLookup1( genome, kmers_plus );
     vec< kmer<K> > kmers( kmers_plus.size( ) );
     for ( int i = 0; i < kmers_plus.isize( ); i++ )
          kmers[i] = kmers_plus[i].first;

     vec<Bool> pairs_to_delete( pairs.nPairs( ), False );
     vec<Bool> reads_to_delete( bases.size( ), False );
     #pragma omp parallel for
     for ( int pid = 0; pid < (int) pairs.nPairs( ); pid++ )
     {    int id1 = pairs.ID1(pid), id2 = pairs.ID2(pid);
          kmer<K> x;
          Bool match = False;
          for ( int j = 0; j <= bases[id1].isize( ) - K; j++ )
          {    x.SetToSubOf( bases[id1], j );
               if ( BinMember( kmers, x ) ) 
               {    pairs_to_delete[pid] = True;
                    reads_to_delete[id1] = reads_to_delete[id2] = True;
                    match = True;
                    break;    }    }
          if (match) continue;
          for ( int j = 0; j <= bases[id2].isize( ) - K; j++ )
          {    x.SetToSubOf( bases[id2], j );
               if ( BinMember( kmers, x ) ) 
               {    pairs_to_delete[pid] = True;
                    reads_to_delete[id1] = reads_to_delete[id2] = True;
                    break;    }    }    }

     // Copied from below.

     vec<int> to_new( bases.size( ) );
     int new_id = 0;
     for ( int id = 0; id < (int) bases.size( ); id++ )
     {    to_new[id] = new_id;
          if ( !reads_to_delete[id] ) new_id++;    }
     pairs.removePairs(pairs_to_delete);
     for ( int pid = 0; pid < (int) pairs.nPairs( ); pid++ )
          pairs.SetIDs( pid, to_new[ pairs.ID1(pid) ], to_new[ pairs.ID2(pid) ] );
     bases.EraseIf(reads_to_delete);
     quals.EraseIf(reads_to_delete);
     String head = TMP + "/frag_reads_orig";
     bases.WriteAll( head + ".fastb" );
     quals.WriteAll( head + ".qualb" );
     pairs.Write( head + ".pairs" );
     if ( IsRegularFile( head + ".qltout" ) )
     {    vec<look_align> aligns; 
          LoadLookAligns( head + ".qltout", aligns );
          vec<Bool> to_delete( aligns.size( ), False );
          for ( int id = 0; id < aligns.isize( ); id++ )
               if ( reads_to_delete[ aligns[id].query_id ] ) to_delete[id] = True;
          EraseIf( aligns, to_delete );
          WriteLookAligns( head + ".qltout", aligns );    }    }

void CleanFosmidReads( const String& TMP )
{    String head = TMP + "/frag_reads_orig";
     vecbasevector bases( head + ".fastb" );
     vecqualvector quals( head + ".qualb" );
     PairsManager pairs;
     pairs.Read( head + ".pairs" );

     vecbasevector junk;
     junk.push_back( FosmidVectorForPool2( ) );
     junk.ReadAll( "/wga/dev/references/Escherichia_coli/genome.fastb", True );
     int nj = junk.size( );
     for ( int i = 0; i < nj; i++ )
          junk.push_back( ReverseComplement( junk[i] ) );

     const int K = 20;

     vec< triple<kmer<K>,int,int> > kmers_plus;
     MakeKmerLookup1( junk, kmers_plus );
     vec< kmer<K> > kmers( kmers_plus.size( ) );
     for ( size_t i = 0; i < kmers.size( ); i++ )
          kmers[i] = kmers_plus[i].first;

     vec<Bool> pairs_to_delete( pairs.nPairs( ), False );
     vec<Bool> reads_to_delete( bases.size( ), False );
     for ( int pid = 0; pid < (int) pairs.nPairs( ); pid++ )
     {    vec<int> id(2);
          id[0] = pairs.ID1(pid), id[1] = pairs.ID2(pid);
          vec<Bool> junky(2, False);
          for ( int j = 0; j < 2; j++ )
          {    kmer<K> x;
               x.SetToSubOf( bases[ id[j] ], 0 );
               int p = BinPosition( kmers, x );
               if ( p >= 0 ) junky[j] = True;    }
          if ( junky[0] && junky[1] ) 
          {    pairs_to_delete[pid] = True;
               reads_to_delete[ id[0] ] = reads_to_delete[ id[1] ] = True;    }    }

     const int L = 80;

     vec< triple<kmer<L>,int,int> > kmers_plus2;
     MakeKmerLookup1( junk, kmers_plus2 );
     vec< kmer<L> > kmers2( kmers_plus2.size( ) );
     for ( size_t i = 0; i < kmers2.size( ); i++ )
          kmers2[i] = kmers_plus2[i].first;
     #pragma omp parallel for
     for ( int pid = 0; pid < (int) pairs.nPairs( ); pid++ )
     {    int id1 = pairs.ID1(pid), id2 = pairs.ID2(pid);
          kmer<L> x;
          Bool match = False;
          for ( int j = 0; j <= bases[id1].isize( ) - L; j++ )
          {    x.SetToSubOf( bases[id1], j );
               if ( BinMember( kmers2, x ) ) 
               {    pairs_to_delete[pid] = True;
                    reads_to_delete[id1] = reads_to_delete[id2] = True;
                    match = True;
                    break;    }    }
          if (match) continue;
          for ( int j = 0; j <= bases[id2].isize( ) - L; j++ )
          {    x.SetToSubOf( bases[id2], j );
               if ( BinMember( kmers2, x ) ) 
               {    pairs_to_delete[pid] = True;
                    reads_to_delete[id1] = reads_to_delete[id2] = True;
                    break;    }    }    }

     vec<int> to_new( bases.size( ) );
     int new_id = 0;
     for ( int id = 0; id < (int) bases.size( ); id++ )
     {    to_new[id] = new_id;
          if ( !reads_to_delete[id] ) new_id++;    }
     pairs.removePairs(pairs_to_delete);
     for ( int pid = 0; pid < (int) pairs.nPairs( ); pid++ )
          pairs.SetIDs( pid, to_new[ pairs.ID1(pid) ], to_new[ pairs.ID2(pid) ] );
     bases.EraseIf(reads_to_delete);
     quals.EraseIf(reads_to_delete);
     bases.WriteAll( head + ".fastb" );
     quals.WriteAll( head + ".qualb" );
     pairs.Write( head + ".pairs" );
     if ( IsRegularFile( head + ".qltout" ) )
     {    vec<look_align> aligns; 
          LoadLookAligns( head + ".qltout", aligns );
          vec<Bool> to_delete( aligns.size( ), False );
          for ( int id = 0; id < aligns.isize( ); id++ )
               if ( reads_to_delete[ aligns[id].query_id ] ) to_delete[id] = True;
          EraseIf( aligns, to_delete );
          WriteLookAligns( head + ".qltout", aligns );    }    }

void RemoveFosmidVector( SupportedHyperBasevector& shb,
     const long_logging_control& log_control, const long_logging& logc,
     const Bool delrc )
{
     double clock = WallClockTime( );
     cout << Date( ) << ": removing Fosmid vector" << endl;

     // Fosmid vector boundaries for hpool2 and hpool3.
     // (tried to put this above in anonymous namespace but something weird
     // happened - sequences were corrupted)

     basevector v1(  "ACACGACGCTCTTCCGATCTAGTTGCTT" ); // genome follows this
     basevector v2(  "GACGTGTGCTCTTCCGATCTAGTTGCTT" ); // genome follows this

     // Load nonstandard junction sites for hpool2 and hpool3.  Junction sites have 
     // the form left|right, where right is the start of the genomic sequence.  
     // Currently left and right are required to each be 20 bases long.

     vec<String> nsj, nsjrc; // nonstandard junctions and their reverse complements
     vec<String> regions;
     vec< vec< pair<String,String> > > junctions, breaks, edits; 
     ParseFosmidPoolMetainfo( regions, junctions, breaks, edits );
     const int ntag = 20;
     for ( int i = 0; i < junctions.isize( ); i++ )
     for ( int j = 0; j < junctions[i].isize( ); j++ )
     {    ForceAssertEq( junctions[i][j].first.isize( ), ntag );
          ForceAssertEq( junctions[i][j].second.isize( ), ntag );
          String jun = junctions[i][j].first + junctions[i][j].second, junrc;
          StringReverseComplement( jun, junrc );
          nsj.push_back(jun), nsjrc.push_back(junrc);    }
     vec< pair<String,String> > all_breaks;
     for ( int i = 0; i < breaks.isize( ); i++ )
     for ( int j = 0; j < breaks[i].isize( ); j++ )
          all_breaks.push_back( breaks[i][j] );

     // Cut at boundaries.

     String s1 = v1.ToString( ), s2 = v2.ToString( ); 
     String s1rc, s2rc;
     StringReverseComplement( s1, s1rc );
     StringReverseComplement( s2, s2rc );
     vec<int> to_left, to_right;
     shb.ToLeft(to_left), shb.ToRight(to_right);
     int ne = shb.EdgeObjectCount( );
     vec<Bool> left_trimmed( ne, False ), right_trimmed( ne, False );
     vec<Bool> trimmed_to_nothing( ne, False );
     vec<int> edels;
     for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
     {    String E = shb.EdgeObject(e).ToString( );
          int left_trim = 0, right_trim = 0;

          for ( int j = 0; j < E.isize( ); j++ )
          {    if ( E.Contains( s1, j ) )
               {    cout << "found left fw at " << e << "." << j << "-"
                         << j + s1.isize( ) << endl;    
                    left_trim = Max( left_trim, j + s1.isize( ) );    }
               if ( E.Contains( s2, j ) )
               {    cout << "found right fw at " << e << "." << j << "-"
                         << j + s2.isize( ) << endl;
                    left_trim = Max( left_trim, j + s2.isize( ) );    }
               if ( E.Contains( s1rc, j ) )
               {    cout << "found left rc at " << e << "." << j << "-"
                         << j + s1.isize( ) << endl;
                    right_trim = Max( right_trim, E.isize( ) - j );    }
               if ( E.Contains( s2rc, j ) )
               {    cout << "found right rc at " << e << "." << j << "-"
                         << j + s2.isize( ) << endl;    
                    right_trim = Max( right_trim, E.isize( ) - j );    }    }

          for ( int j = 0; j < E.isize( ); j++ )
          for ( int l = 0; l < nsj.isize( ); l++ )
          {    if ( E.Contains( nsj[l], j ) )
               {    cout << "found nonstandard junction " << l << " fw at "
                         << e << "." << j + ntag << endl;
                    left_trim = Max( left_trim, j + ntag );    }
               if ( E.Contains( nsjrc[l], j ) )
               {    cout << "found nonstandard junction " << l << " rc at "
                         << e << "." << E.isize( ) - j - ntag << endl;
                    right_trim = Max( right_trim, E.isize( ) - j - ntag );    }    }

          // Process breakpoints.  Currently we just keep the longer edge, not
          // really right.

          vec<int> breakpoints;
          for ( int p = 0; p < E.isize( ); p++ )
          {    for ( int j = 0; j < all_breaks.isize( ); j++ )
               {    String s1 = all_breaks[j].first, s2 = all_breaks[j].second;
                    String s1rc, s2rc;
                    StringReverseComplement( s1, s1rc );
                    StringReverseComplement( s2, s2rc );
                    if ( E.Contains( s1 + s2, p ) )
                    {    breakpoints.push_back( p + s1.isize( ) );
                         cout << "found breakpoint at " << e << "." << p 
                              << endl;    }
                    if ( E.Contains( s2rc + s1rc, p ) )
                    {    breakpoints.push_back( p + s2rc.isize( ) );
                         cout << "found breakpoint at " << e << "." << p 
                              << endl;    }    }    }
          UniqueSort(breakpoints);
          if ( breakpoints.solo( ) )
          {    int p = breakpoints[0];
               String left = E.substr( 0, p ), right = E.substr( p, E.isize( ) - p );
               if ( left.size( ) > right.size( ) )
                    right_trim = Max( right_trim, E.isize( ) - p );
               if ( right.isize( ) > left.isize( ) )
                    left_trim = Max( left_trim, p );    }

          int new_size = shb.EdgeObject(e).isize( ) - left_trim - right_trim;
          if ( new_size >= shb.K( ) )
          {    shb.EdgeObjectMutable(e).SetToSubOf( 
                    shb.EdgeObject(e), left_trim, new_size );    }
          else 
          {    trimmed_to_nothing[e] = True;
               edels.push_back(e);    }
          int v = to_left[e], w = to_right[e];
          if ( left_trim > 0 && shb.To(v).size( ) + shb.From(v).size( ) > 1 )
          {    shb.AddVertices(1);
               shb.GiveEdgeNewFromVx( e, v, shb.N( ) - 1 );
               left_trimmed[e] = True;    }
          if ( right_trim > 0 && shb.To(w).size( ) + shb.From(w).size( ) > 1 )
          {    shb.AddVertices(1);
               shb.GiveEdgeNewToVx( e, w, shb.N( ) - 1 );
               right_trimmed[e] = True;    }    }
     vec<Bool> to_delete( shb.NPaths( ), False );
     for ( int i = 0; i < shb.NPaths( ); i++ )
     {    for ( int j = 0; j < shb.Path(i).isize( ); j++ )
          {    int e = shb.Path(i)[j];
               if ( trimmed_to_nothing[e] || ( left_trimmed[e] && j > 0 ) )
                    to_delete[i] = True;
               if ( trimmed_to_nothing[e] 
                    && ( right_trimmed[e] && j < shb.Path(i).isize( ) - 1 ) )
               {    to_delete[i] = True;    }    }    }
     EraseIf( shb.PathsMutable( ), to_delete );
     // Temporary if.
     if ( shb.WeightsFwOrigin( ).size( ) == shb.WeightsFw( ).size( ) )
     {    EraseIf( shb.WeightsFwOriginMutable( ), to_delete );
          EraseIf( shb.WeightsRcOriginMutable( ), to_delete );    }
     EraseIf( shb.WeightsFwMutable( ), to_delete );
     EraseIf( shb.WeightsRcMutable( ), to_delete );

     shb.DeleteEdges(edels);
     shb.RemoveUnneededVertices( );
     shb.RemoveDeadEdgeObjects( );
     shb.RemoveEdgelessVertices( );

     // Carry out edits.  Note that this could break the structure.  Also terrible
     // since we're ignoring the id of the Fosmid.

     for ( int i = 0; i < edits.isize( ); i++ )
     for ( int l = 0; l < edits[i].isize( ); l++ )
     {    String x = edits[i][l].first, y = edits[i][l].second, xrc, yrc;
          StringReverseComplement( x, xrc );
          StringReverseComplement( y, yrc );
          for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
          {    String E = shb.EdgeObject(e).ToString( );
               for ( int j = 0; j < shb.EdgeObject(e).isize( ); j++ )
               {    if ( E.Contains( x, j ) )
                    {    E.ReplaceBy( x, y );
                         shb.EdgeObjectMutable(e) = basevector(E);    }
                    if ( E.Contains( xrc, j ) )
                    {    E.ReplaceBy( xrc, yrc );
                         shb.EdgeObjectMutable(e) 
                              = basevector(E);    }    }    }    }

     // Check for non-adjacent edges in paths.  This most be compensating for
     // a bug in the above code.

     shb.ToRight(to_right);
     vec<Bool> pdel( shb.NPaths( ), False );
     for ( int i = 0; i < shb.NPaths( ); i++ )
     {    for ( int j = 0; j < shb.Path(i).isize( ) - 1; j++ )
          {    int e1 = shb.Path(i)[j], e2 = shb.Path(i)[j+1];
               Bool adj = False;
               int v = to_right[e1];
               for ( int l = 0; l < shb.From(v).isize( ); l++ )
                    if ( shb.EdgeObjectIndexByIndexFrom( v, l ) == e2 ) adj = True;
               if ( !adj ) pdel[i] = True;    }    }
     EraseIf( shb.PathsMutable( ), pdel );
     // Temporary if.
     if ( shb.WeightsFwOrigin( ).size( ) == shb.WeightsFw( ).size( ) )
     {    EraseIf( shb.WeightsFwOriginMutable( ), pdel );
          EraseIf( shb.WeightsRcOriginMutable( ), pdel );    }
     EraseIf( shb.WeightsFwMutable( ), pdel );
     EraseIf( shb.WeightsRcMutable( ), pdel );

     REPORT_TIME( clock, "used removing Fosmid vector" );
     if (delrc) shb.DeleteReverseComplementComponents(logc);     }

void FlagVector( const VecEFasta& creads, vec<Bool>& to_delete,
     const long_logging& logc )
{
     // Identify cloning site.

     double clock = WallClockTime( );
     const int L = 40;
     vec<basevector> heads;
     heads.push_back( vleft, vright );
     vec<basevector> hs;
     for ( int j = 0; j < heads.isize( ); j++ )
     {    vec<basevector> posts, postsx;
          #pragma omp parallel for
          for ( size_t id = 0; id < creads.size( ); id++ )
          {    const efasta& e = creads[id];
               vec<basevector> x;
               e.ExpandTo(x);
               for ( int l = 0; l < x.isize( ); l++ )
               {    const basevector& r = x[l];
                    String rs = r.ToString( );
                    for ( int pos = 0; pos <= r.isize( ) - heads[j].isize( ) - L; 
                         pos++ )
                    {    if ( rs.Contains( heads[j].ToString( ), pos ) )
                         {    
                              #pragma omp critical
                              {    posts.push( r, pos + heads[j].isize( ), 
                                        L );    }    }    }    }    }
          UniqueSort(posts);
          vec<int> freq;
          for ( int i = 0; i < posts.isize( ); i++ )
          {    int j = posts.NextDiff(i);
               postsx.push_back( posts[i] );
               freq.push_back( j - i );
               i = j - 1;    }
          ReverseSortSync( freq, postsx );
          if ( posts.nonempty( ) ) hs.push_back( Cat( heads[j], postsx[0] ) );    }

     // Identify vector.

     basevector fvec = FosmidVector( );
     vecbasevector v;
     v.push_back(fvec);
     fvec.ReverseComplement( );
     v.push_back(fvec);
     for ( int j = 0; j < hs.isize( ); j++ )
     {    v.push_back( hs[j] );
          hs[j].ReverseComplement( );
          v.push_back( hs[j] );    }
     vec< triple<kmer<L>,int,int> > kmers_plus;
     MakeKmerLookup1( v, kmers_plus );
     vec< kmer<L> > kmers( kmers_plus.size( ) );
     for ( size_t i = 0; i < kmers.size( ); i++ )
          kmers[i] = kmers_plus[i].first;
     vec<Bool> evil( creads.size( ), False );
     #pragma omp parallel for
     for ( int id = 0; id < (int) creads.size( ); id++ )
     {    const efasta& e = creads[id];
          vec<basevector> x;
          e.ExpandTo(x);
          for ( int l = 0; l < x.isize( ); l++ )
          {    const basevector& r = x[l];
               kmer<L> x;
               for ( int pos = 0; pos <= r.isize( ) - L; pos++ )
               {    x.SetToSubOf( r, pos );
                    if ( BinMember( kmers, x ) )
                    {    to_delete[id] = True;
                         break;    }    }    }    }
     REPORT_TIME( clock, "used flagging vector" );    }
