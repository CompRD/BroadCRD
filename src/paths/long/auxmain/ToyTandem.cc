///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// note results changed now should rerun example assemblies

// ToyTandem.  Play with tandem repeat assembly.
//
// Test case 1.
// r45873:LongProto SAMPLE=human READS=#picard TMP=tmp.ttt1 OUT_INT_HEAD=ttt1
//     DATA_SPEC=HUMAN_CONTROLS=52 HEURISTICS=ORIENT_TO_REFERENCE=True
// ToyTandem IN_SHBV=ttt1.final.shbv TMP=tmp.ttt1 BOUNDS=129,127 MAX_DUP=2
// truth = 129,126,128,126,121,21,20,130,23,20,22,121,127
//
// Test case 2.
// r45873:LongProto SAMPLE=rhody X=5 READS=#picard TMP=tmp.ttt2 OUT_INT_HEAD=ttt2
// ToyTandem IN_SHBV=ttt2.final.shbv TMP=tmp.ttt2 BOUNDS=87,73 MAX_DUP=2
// truth = 87,66,86,37,68,39,73
//
// Test case 3.
// r45882:LongProto SAMPLE=human READS=#picard TMP=tmp.ttt3 OUT_INT_HEAD=ttt3
//     DATA_SPEC=HUMAN_CONTROLS=5 HEURISTICS=ORIENT_TO_REFERENCE=True
// ToyTandem IN_SHBV=ttt3.final.shbv TMP=tmp.ttt3 BOUNDS=90,22 MAX_DUP=2
// truth = 90,89,82,89,22
// two concentric loops
//
// Test case 4.
// r45882:LongProto SAMPLE=human READS=#picard TMP=tmp.ttt4 OUT_INT_HEAD=ttt4
//     DATA_SPEC=HUMAN_CONTROLS=40 HEURISTICS=ORIENT_TO_REFERENCE=True
// ToyTandem IN_SHBV=ttt4.final.shbv TMP=tmp.ttt4 BOUNDS=31,24 MAX_DUP=2
// true path not present
//
// Test case 5.
// r45882:LongProto SAMPLE=human READS=#picard TMP=tmp.ttt5 OUT_INT_HEAD=ttt5
//     DATA_SPEC=HUMAN_CONTROLS=93 HEURISTICS=ORIENT_TO_REFERENCE=True
// ToyTandem IN_SHBV=ttt5.final.shbv TMP=tmp.ttt5 BOUNDS=83,3 MAX_DUP=2
// true path not present
//
// Test case 6.
// r45873:LongProto SAMPLE=rhody X=5 READS=#picard TMP=tmp.ttt2 OUT_INT_HEAD=ttt2
// ToyTandem IN_SHBV=ttt2.final.shbv TMP=tmp.ttt2 BOUNDS=73,88 MAX_DUP=4
// (runs for a long time)
// truth = 73,16,53,48,35,55,15,59,60,5,0,44,61,64,85,44,61,65,8,55,15,58,3,4,83,0,
//         44,36,52,47,72,14,15,58,3,4,83,0,44,61,63,54,15,58,3,4,83,0,44,36,52,56,88

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "ParseSet.h"
#include "Qualvector.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/EvalByReads.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/SupportedHyperBasevector.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(IN_SHBV, "shbv file for assembly");
     CommandArgument_String_Doc(TMP, "pre-existing tmp dir generated by LongProto");
     CommandArgument_String_Doc(BOUNDS, "left-edge,right-edge, bounding cell");
     CommandArgument_Int_Doc(MAX_DUP, "temporary heuristic - maximum "
          "number of times that edge may be used");
     CommandArgument_Int_OrDefault(VERBOSITY, 1);
     CommandArgument_String_OrDefault_Doc(TRACE_IDS, "", "trace these reads");
     EndCommandArguments;

     // Load data.

     SupportedHyperBasevector shb;
     BinaryReader::readFile( IN_SHBV, &shb );
     vec<int> seq;
     ParseIntSet( BOUNDS, seq, false );
     vecbasevector bases( TMP + "/frag_reads_orig.fastb" );
     vecqualvector quals( TMP + "/frag_reads_orig.qualb" );
     PairsManager pairs;
     pairs.Read( TMP + "/frag_reads_orig.pairs" );

     // Find all paths through the cell.

     vec< vec<int> > paths;
     vec<int> to_left, to_right;
     shb.ToLeft(to_left), shb.ToRight(to_right);
     shb.EdgePaths( to_left[ seq[0] ], to_right[ seq[1] ], 
          paths, MAX_DUP );
     vec<Bool> to_delete( paths.size( ), False );
     for ( int i = 0; i < paths.isize( ); i++ )
     {    const vec<int>& p = paths[i];
          if ( p.size( ) < 2 || p.front( ) != seq[0] || p.back( ) != seq[1] )
          {    to_delete[i] = True;    }    }
     EraseIf( paths, to_delete );
     if ( VERBOSITY >= 2 )
     {    cout << "paths:\n";
          for ( int i = 0; i < paths.isize( ); i++ )
               cout << "[" << i+1 << "]: " << printSeq( paths[i] ) << "\n";    }
     if ( VERBOSITY == 1 ) 
          cout << "There are " << paths.size( ) << " paths." << endl;
     vec<basevector> bpaths( paths.size( ) );
     for ( int i = 0; i < paths.isize( ); i++ )
          bpaths[i] = shb.Cat( paths[i] );
     vec<int> middle;
     for ( int i = 0; i < paths.isize( ); i++ )
     for ( int j = 1; j < paths[i].isize( ) - 1; j++ )
          middle.push_back( paths[i][j] );
     UniqueSort(middle);

     // Find the pairs that map to the locus.

     vec<int> trace_ids;
     ParseIntSet( TRACE_IDS, trace_ids );
     const int L = 12;
     HyperBasevector hb_fw(shb), hb_rc(shb);
     hb_rc.Reverse( );
     vec<int> to_right_fw, to_right_rc;
     hb_fw.ToRight(to_right_fw), hb_rc.ToRight(to_right_rc);
     vecbasevector x_fw, x_rc;
     for ( int i = 0; i < hb_fw.EdgeObjectCount( ); i++ )
          x_fw.push_back( hb_fw.EdgeObject(i) );
     for ( int i = 0; i < hb_rc.EdgeObjectCount( ); i++ )
          x_rc.push_back( hb_rc.EdgeObject(i) );
     VecIntPairVec locs_fw, locs_rc;
     CreateGlocs( x_fw, L, locs_fw );
     CreateGlocs( x_rc, L, locs_rc );
     vec< vec<read_place> > PLACES( bases.size( ) );

     Bool filter = True;
     const int prox = 20;

     // Bool filter = False;
     // const int prox = 0;

     #pragma omp parallel for
     for ( int64_t id = 0; id < (int64_t) bases.size( ); id++ )
     {    int n = KmerId( bases[id], L, 0 );
          const int infinity = 1000000000;
          int qual_sum = infinity;
          const int min_qual = 1;
          FindPlaces( bases[id], quals[id], n, hb_fw, hb_rc, to_right_fw,
               to_right_rc, locs_fw, locs_rc, PLACES[id], 
               qual_sum, min_qual, prox * 1000 );

          if ( BinMember( trace_ids, id ) )
          {
               #pragma omp critical
               {    cout << "\noriginal placements of read " << id << "\n";
                    for ( int j = 0; j < PLACES[id].isize( ); j++ )
                    {    cout << "[" << j+1 << "] " << PLACES[id][j] 
                              << endl;    }    }    }

          if ( !filter ) continue;

          // Reduce to a single unique placement.  Note that this does not
          // correctly set the score.

          if ( PLACES[id].size( ) > 1 )
          {    
               // Raise quality scores, then rescore alignments and filter by score.

               const int flank = 10;
               const basevector& b = bases[id];
               qualvector q = quals[id];
               for ( int j = 0; j < PLACES[id].isize( ); j++ )
               {    const read_place& p = PLACES[id][j];
                    const HyperBasevector& hb = ( p.Fw( ) ? hb_fw : hb_rc );
                    int ei = 0, pos = p.P( ), matches = 0;
                    for ( int l = 0; l < b.isize( ); l++ )
                    {    if ( b[l] != hb.EdgeObject( p.E(ei) )[pos] ) matches = 0;
                         else matches++;
                         if ( matches >= 2*flank+1 )
                              q[l-flank] = Max( (int) q[l-flank], prox );
                         pos++;
                         if ( pos == hb.EdgeObject( p.E(ei) ).isize( ) )
                         {    ei++;
                              if ( ei == p.N( ) ) break;
                              pos = hb.K( ) - 1;    }    }    }
               int minq = 1000000000;
               for ( int j = 0; j < PLACES[id].isize( ); j++ )
               {    read_place& p = PLACES[id][j];
                    p.ComputeQsum( b, q, ( p.Fw( ) ? hb_fw : hb_rc ), min_qual );
                    minq = Min( minq, p.Qsum( ) );    }
               vec<Bool> to_delete( PLACES[id].size( ), False );
               for ( int j = 0; j < PLACES[id].isize( ); j++ )
               {    if ( PLACES[id][j].Qsum( ) >= minq + prox * 1000 ) 
                         to_delete[j] = True;    }
               EraseIf( PLACES[id], to_delete );
               if ( PLACES[id].solo( ) ) continue;

               if ( BinMember( trace_ids, id ) )
               {
                    #pragma omp critical
                    {    cout << "\nprefiltered placements of read " << id << "\n";
                         for ( int j = 0; j < PLACES[id].isize( ); j++ )
                         {    cout << "[" << j+1 << "] " << PLACES[id][j] 
                                   << endl;    }    }    }

               Bool fw = False, rc = False;
               vec<int> p;
               for ( int j = 0; j < PLACES[id].isize( ); j++ )
               {    if ( PLACES[id][j].Fw( ) ) fw = True;
                    else rc = True;
                    p.push_back( PLACES[id][j].P( ) );    }
               UniqueSort(p);
               if ( ( fw && rc ) || p.size( ) > 1 )
               {    PLACES[id].clear( );
                    continue;    }
               int k;
               for ( k = 0; ; k++ )
               {    Bool done = False;
                    for ( int j = 0; j < PLACES[id].isize( ); j++ )
                    {    if ( PLACES[id].isize( ) == k )
                         {    done = True;
                              break;    }
                         for ( int l = 1; l < PLACES[id].isize( ); l++ )
                         {    if ( PLACES[id][l].E(k) != PLACES[id][0].E(k) )
                              {    done = True;
                                   break;    }    }    }
                    if (done) break;    }
               if ( k == 0 )
               {    PLACES[id].clear( );
                    continue;    }
               PLACES[id].resize(1);
               PLACES[id][0].EMutable( ).resize(k);    }    }
     // cout << Date( ) << ": done" << endl;
     // cout << "\n";

     vec<int64_t> pids;
     for ( int64_t pid = 0; pid < (int64_t) pairs.nPairs( ); pid++ )
     {    int64_t id1 = pairs.ID1(pid), id2 = pairs.ID2(pid);
          Bool good = False;
          for ( int pass = 0; pass < 2; pass++ )
          {    int64_t id = ( pass == 0 ? id1 : id2 );
               for ( int j = 0; j < PLACES[id].isize( ); j++ )
               for ( int k = 0; k < PLACES[id][j].N( ); k++ )
                    if ( BinMember( middle, PLACES[id][j].E(k) ) ) good = True;    }
          if (good) pids.push_back(pid);    }
     cout << "There are " << pids.size( ) << " pairs incident upon the cell." << endl;

     for ( int j = 0; j < paths.isize( ); j++ )
     {    vec<String> cov;
          if ( VERBOSITY >= 2 )
               cout << "\npath " << j+1 << " = " << printSeq( paths[j] ) << endl;
          for ( int i = 0; i < pids.isize( ); i++ )
          {    int64_t id1 = pairs.ID1( pids[i] ), id2 = pairs.ID2( pids[i] );
               for ( int l1 = 0; l1 < PLACES[id1].isize( ); l1++ )
               for ( int l2 = 0; l2 < PLACES[id2].isize( ); l2++ )
               {    read_place p1 = PLACES[id1][l1], p2 = PLACES[id2][l2];
                    int64_t jd1(id1), jd2(id2);
                    if ( !p1.Fw( ) ) 
                    {    swap( p1, p2 );
                         swap( jd1, jd2 );    }
                    if ( !p1.Fw( ) || p2.Fw( ) ) continue;
                    p2.Reverse( bases[jd2], shb );
                    for ( int u1 = 0; u1 < paths[j].isize( ); u1++ )
                    {    if ( paths[j].Contains( p1.E( ), u1 ) )
                         {    for ( int u2 = u1; u2 < paths[j].isize( ); u2++ )
                              {    if ( paths[j].Contains( p2.E( ), u2 ) )
                                   {    String c( paths[j].size( ), ' ' );
                                        for ( int l = 0; l < p1.N( ); l++ )
                                             c[ u1 + l ] = '-';
                                        for ( int l = 0; l < p2.N( ); l++ )
                                             c[ u2 + l ] = '-';
                                        cov.push_back(c);
                                        if ( VERBOSITY >= 2 )
                                        {
                                        cout << jd1 << "/" << jd2 << " --> " 
                                             << p1.E(0) << "." << p1.P( );
                                        for ( int j = 1; j < p1.N( ); j++ )
                                             cout << "," << p1.E(j);
                                        cout << " [" << setiosflags(ios::fixed) 
                                             << setprecision(1)
                                             << p1.Qsum( ) / 1000.0 
                                             << resetiosflags(ios::fixed) << "]";
                                        cout << " :: ";
                                        cout << p2.E(0) << "." << p2.P( );
                                        for ( int j = 1; j < p2.N( ); j++ )
                                             cout << "," << p2.E(j);
                                        cout << " [" << setiosflags(ios::fixed) 
                                             << setprecision(1)
                                             << p2.Qsum( ) / 1000.0 
                                             << resetiosflags(ios::fixed) << "]"
                                             << "\n";    
                                        }
                                        }    }    }    }    }    }

          // Display coverage.

          Sort(cov);
          vec<Bool> gapped( cov.size( ), False );
          for ( int l = 0; l < cov.isize( ); l++ )
               if ( cov[l].After( "-" ).Contains( " -" ) ) gapped[l] = True;
          int mc = 1000000000, mcc = 1000000000;
          for ( int u = 0; u < paths[j].isize( ); u++ )
          {    int c = 0, cc = 0;;
               for ( int l = 0; l < cov.isize( ); l++ )
               {    if ( cov[l][u] == '-' ) c++;
                    if ( !gapped[l] && cov[l][u] == '-' ) cc++;    }
               mc = Min( mc, c );
               mcc = Min( mcc, cc );    }
          if ( mc <= 1 && VERBOSITY <= 1 ) continue;
          if ( VERBOSITY == 1 )
               cout << "\npath " << j+1 << " = " << printSeq( paths[j] ) << endl;
          cout << "\n";
          for ( int l = 0; l < cov.isize( ); l++ )
          {    int m;
               for ( m = l + 1; m < cov.isize( ); m++ )
                    if ( cov[m] != cov[l] ) break;
               cout << cov[l] << "  " << m-l << "\n";    
               l = m - 1;    }
          cout << "\nminimum coverage = " << mc << endl;
          cout << "minimum contiguous coverage = " << mcc << endl;
          cout << "ratio = " << setiosflags(ios::fixed) << setprecision(1)
               << double(mc)/double(mcc) << resetiosflags(ios::fixed) 
               << endl;    }    }
