///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Experimental code for recovering kmers that sequence poorly.

// Note experimental versions Snorgle1.cc and Snorgle2.cc need to be kept in sync 
// with phases 1 and 2 of this.

// Find read-read alignments for reads that might be in gaps.

// PHASE 1
//
// Design of finding read-read alignments:
// 1. Use seed kmers occurring no more than max_mult times.
// 2. For a given perfect match between two reads, with the first one forward,
//    build an alignment only from the leftmost kmer in the perfect match.
// 3. However, condition #1 still applies to all matching kmers.
// Getting these to work compatibly is tricky.
//
// Notes.
// 1. Does not handle palindrome seeds correctly.  This is because the kmer record
//    could (in principle) be twice in the stack in two forms, but isn't.  Can be 
//    remedied by adding in flipped alignments.  Most of these palindromes are AT
//    homopolymers.
// 2. Seems to duplicate calculations two-fold.
// 3. Main loop ( @ third Process stack ) doesn't really have to be n^2.
// 4. Wondering if homopolymer/dinuke-based stacks greatly increase run time.
// 5. Could reduce memory by running in 4 passes.
// 6. SortSync convoluted.

// GapToy X=fos2 PAD=-10K INSTANCE=s

// GapToy X=fos2 PAD=1M INSTANCE=1

// GapToy X=fos1 PAD=10M INSTANCE=2 (17 min crd8)
// nfriends = 9262671
// 22.3 seconds (crd8) (outdated)

// GapToy X=fos3 PAD=20M INSTANCE=3 (28 min crd28)
// nfriends = 14165639
// 36.7 seconds (crd8)
// scales to 0.82 hours

// GapToy X=fos3 PAD=40M INSTANCE=4 (45 min crd28)
// nfriends = 24393063
// 67 seconds (crd8)
// 37 seconds (crd28)

// GapToy X=fos3 PAD=80M INSTANCE=5 (1.77 hours crd28)
// nfriends = 135659521
// 7.0 minutes (crd28)
// peak mem = 83.5 GB

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "MainTools.h"
#include "Qualvector.h"
#include "VecUtilities.h"
#include "feudal/PQVec.h"
#include "kmers/KmerRecord.h"
#include "kmers/MakeLookup.h"
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/FriendAligns.h"
#include "paths/long/LoadCorrectCore.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/ReadPath.h"
#include "paths/long/ReadStack.h"
#include "paths/long/large/explore/SnorgleTools.h"

namespace { // open anonymous namespace

Bool HQDiff( const int64_t id1, const int64_t id2, const int offset, const Bool fw2,
     const vecbasevector& B, const vecqualvector& Q )
{
     const int qtop = 30;
     int low = Max( 0, offset );
     if (fw2)
     {    for ( int p1 = low; p1 < B[id1].isize( ); p1++ )
          {    int p2 = p1 - offset;
               if ( p2 >= B[id2].isize( ) ) break;
               if ( B[id1][p1] != B[id2][p2] )
               {    if ( Q[id1][p1] >= qtop && Q[id2][p2] >= qtop ) 
                         return True;    }    }    }
     else
     {    int n2 = B[id2].size( );
          for ( int p1 = low; p1 < B[id1].isize( ); p1++ )
          {    int p2 = p1 - offset;
               if ( p2 >= B[id2].isize( ) ) break;
               if ( B[id1][p1] != 3 - B[id2][ n2 - p2 - 1 ] )
               {    if ( Q[id1][p1] >= qtop && Q[id2][ n2 - p2 - 1 ] >= qtop ) 
                         return True;    }    }    }
     return False;    }

template<int K>
void ProcessRight( vec< triple<kmer<K>,int,int> >& right, const vecbasevector& B,
     const vecqualvector& Q, const int max_mult, 
     vec< vec< quad<int64_t,int64_t,int,Bool> > >& friendsi, 
     const int bi, const vec<int>& ids, const int TRACE_PID )
{    
     vec<int> orig_pos( right.size( ) );
     for ( int64_t j = 0; j < right.jsize( ); j++ )
          orig_pos[j] = right[j].third;

     for ( int64_t m = 0; m < right.jsize( ); m++ )
     {    restart:

          int64_t n;
          for ( n = m + 1; n < right.jsize( ); n++ )
               if ( right[n].first != right[m].first ) break;

          if ( n - m == 1 ) continue;
          if ( n - m > max_mult )
          {    for ( int64_t k = m; k < n; k++ )
               {    int id = right[k].second;
                    int& pos = right[k].third;
                    if ( pos >= 0 )
                    {    pos++;
                         if ( pos + K > B[id].isize( ) )
                         {    if ( k > m ) 
                              {    swap( right[k], right[m] );
                                   swap( orig_pos[k], orig_pos[m] );    }
                              m++;    }
                         else right[k].first.SetToSubOf( B[id], pos );    }
                    else
                    {    pos++;
                         if ( pos == 0 )
                         {    if ( k > m )
                              {    swap( right[k], right[m] );
                                   swap( orig_pos[k], orig_pos[m] );    }
                              m++;    }
                         else
                         {    right[k].first.SetToSubOf( B[id], -pos-1 );
                              right[k].first.ReverseComplement( );    }    }    }

               vec< triple<kmer<K>,int,int> > r(n-m);
               vec<int> rp(n-m);
               for ( int64_t l = m; l < n; l++ )
               {    r[l-m] = right[l];
                    rp[l-m] = orig_pos[l];    }
               SortSync( r, rp );
               for ( int64_t l = m; l < n; l++ )
               {    right[l] = r[l-m];
                    orig_pos[l] = rp[l-m];    }
               
               // SortSync( right.begin( ) + m, right.begin( ) + n, orig_pos );

               goto restart;    }

          // Process stack.

          for ( int64_t k1 = m; k1 < n; k1++ )
          {    int pos1 = right[k1].third; 
               if ( pos1 < 0 ) continue;
               for ( int64_t k2 = m; k2 < n; k2++ )
               {    if ( k1 == k2 ) continue;
                    int id1 = right[k1].second, id2 = right[k2].second;
                    if ( id1 == id2 ) continue;
                    int pos2 = right[k2].third;
                    const basevector &b1 = B[id1], &b2 = B[id2];

                    // Make sure we were at left.

                    int opos1 = orig_pos[k1], opos2 = orig_pos[k2];
                    if ( opos1 >= 0 && opos2 >= 0 )
                    {    if ( opos1 > 0 && opos2 > 0 && b1[opos1-1] == b2[opos2-1] )
                              continue;    }
                    else if ( opos1 < 0 && opos2 < 0 )
                    {    if ( -opos1-1 > 0 && -opos2-1 > 0
                              && b1[-opos1-2] == b2[-opos2-2] ) 
                         {    continue;    }    }
                    else if ( opos1 >= 0 && opos2 < 0 )
                    {    if ( opos1 > 0 && -opos2+K <= b2.isize( )
                              && b1[opos1-1] == 3 - b2[-opos2+K-1] )
                         {    continue;    }    }
                    else if ( opos1 < 0 && opos2 >= 0 )
                    {    if ( -opos1-1 > 0 && opos2 + K < b2.isize( )
                              && 3 - b1[-opos1-2] == b2[opos2+K] )
                         {    continue;    }    }

                    // Trace.
     
                    if ( TRACE_PID >= 0 && TRACE_PID == ids[id1] / 2 )
                    {
                         #pragma omp critical
                         {    cout << "see overlap of " << ids[id1] << " with "
                                   << ids[id2] << ", pos1 = " << pos1
                                   << ", pos2 = " << pos2 << endl;    }    }
                    if ( TRACE_PID >= 0 && TRACE_PID == ids[id2] / 2 )
                    {
                         #pragma omp critical
                         {    cout << "see overlap of " << ids[id2] << " with "
                                   << ids[id1] << ", pos1 = " << pos1
                                   << ", pos2 = " << pos2 << endl;    }    }
     
                    // Form alignment.
     
                    Bool fw2 = !( pos2 < 0 );
                    int xpos1 = pos1, xpos2 = pos2;
                    if ( pos2 < 0 ) xpos2 = -xpos2-1;
                    if ( !fw2 ) xpos2 = B[id2].isize( ) - xpos2 - K;
                    int offset = xpos1 - xpos2;

                    // Kill if hq diff, then save.

                    if ( HQDiff( id1, id2, offset, fw2, B, Q ) ) continue;
                    friendsi[bi].push( id1, id2, offset, fw2 );    }    }
          m = n - 1;    }    }

template<int K>
void ProcessLeft( vec< triple<kmer<K>,int,int> >& left, const vecbasevector& B,
     const vecqualvector& Q, const int max_mult, 
     vec< vec< quad<int64_t,int64_t,int,Bool> > >& friendsi, 
     const int bi, const vec<int>& ids, const int TRACE_PID )
{    
     vec<int> orig_pos( left.size( ) );
     for ( int64_t j = 0; j < left.jsize( ); j++ )
          orig_pos[j] = left[j].third;

     for ( int64_t m = 0; m < left.jsize( ); m++ )
     {    restart:

          int64_t n;
          for ( n = m + 1; n < left.jsize( ); n++ )
               if ( left[n].first != left[m].first ) break;
          if ( n - m == 1 ) continue;
          if ( n - m > max_mult )
          {    for ( int64_t k = m; k < n; k++ )
               {    int id = left[k].second;
                    int& pos = left[k].third;
                    if ( pos >= 0 )
                    {    if ( pos == 0 )
                         {    if ( k > m ) 
                              {    swap( left[k], left[m] );
                                   swap( orig_pos[k], orig_pos[m] );    }
                              m++;    }
                         else
                         {    pos--;
                              left[k].first.SetToSubOf( B[id], pos );    }    }
                    else
                    {    pos--;
                         if ( -pos-1 + K > B[id].isize( ) )
                         {    if ( k > m )
                              {    swap( left[k], left[m] );
                                   swap( orig_pos[k], orig_pos[m] );    }
                              m++;    }
                         else
                         {    left[k].first.SetToSubOf( B[id], -pos-1 );
                              left[k].first.ReverseComplement( );    }    }    }

               vec< triple<kmer<K>,int,int> > r(n-m);
               vec<int> rp(n-m);
               for ( int64_t l = m; l < n; l++ )
               {    r[l-m] = left[l];
                    rp[l-m] = orig_pos[l];    }
               SortSync( r, rp );
               for ( int64_t l = m; l < n; l++ )
               {    left[l] = r[l-m];
                    orig_pos[l] = rp[l-m];    }
               
               // SortSync( left.begin( ) + m, left.begin( ) + n, orig_pos );

               goto restart;    }

          // Process stack.

          for ( int64_t k1 = m; k1 < n; k1++ )
          {    int pos1 = left[k1].third; 
               if ( pos1 >= 0 ) continue;
               for ( int64_t k2 = m; k2 < n; k2++ )
               {    if ( k1 == k2 ) continue;
                    int id1 = left[k1].second, id2 = left[k2].second;
                    if ( id1 == id2 ) continue;
                    int pos2 = left[k2].third;
                    const basevector &b1 = B[id1], &b2 = B[id2];

                    // Make sure we were at left.

                    int opos1 = orig_pos[k1], opos2 = orig_pos[k2];
                    if ( opos1 >= 0 && opos2 >= 0 )
                    {    if ( opos1 > 0 && opos2 > 0 && b1[opos1-1] == b2[opos2-1] )
                              continue;    }
                    else if ( opos1 < 0 && opos2 < 0 )
                    {    if ( -opos1-1 > 0 && -opos2-1 > 0
                              && b1[-opos1-2] == b2[-opos2-2] ) 
                         {    continue;    }    }
                    else if ( opos1 >= 0 && opos2 < 0 )
                    {    if ( opos1 > 0 && -opos2+K <= b2.isize( )
                              && b1[opos1-1] == 3 - b2[-opos2+K-1] )
                         {    continue;    }    }
                    else if ( opos1 < 0 && opos2 >= 0 )
                    {    if ( -opos1-1 > 0 && opos2 + K < b2.isize( )
                              && 3 - b1[-opos1-2] == b2[opos2+K] )
                         {    continue;    }    }
     
                    // Trace.
     
                    if ( TRACE_PID >= 0 && TRACE_PID == ids[id1] / 2 )
                    {
                         #pragma omp critical
                         {    cout << "see overlap of " << ids[id1] << " with "
                                   << ids[id2] << ", pos1 = " << pos1
                                   << ", pos2 = " << pos2 << endl;    }    }
                    if ( TRACE_PID >= 0 && TRACE_PID == ids[id2] / 2 )
                    {
                         #pragma omp critical
                         {    cout << "see overlap of " << ids[id2] << " with "
                                   << ids[id1] << ", pos1 = " << pos1
                                   << ", pos2 = " << pos2 << endl;    }    }
     
                    // Form alignment.
     
                    Bool fw2 = ( pos2 < 0 );
                    int xpos1 = -pos1-1;
                    int xpos2 = pos2;
                    if ( pos2 < 0 ) xpos2 = -pos2-1;
                    if ( !fw2 ) xpos2 = B[id2].isize( ) - xpos2 - K;
                    int offset = xpos1 - xpos2;

                    // Kill if hq diff, then save.

                    if ( HQDiff( id1, id2, offset, fw2, B, Q ) ) continue;
                    friendsi[bi].push( id1, id2, offset, fw2 );    }    }
          m = n - 1;    }    }

} // close anonymous namespace

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(PHASE, "0 or 1 or 2 or 3 or 4 or list of these");
     CommandArgument_Int_OrDefault_Doc(TRACE_PID, -1,
          "for phase 2, analyze just this PID; also follow in phase 1");
     CommandArgument_Int_OrDefault_Doc(TRACE_RID, -1,
          "sets TRACE_PID = TRACE_RID/2");
     CommandArgument_String_OrDefault_Doc(INSTANCE, "51400.newchem",
          "s or 1 or 2 or 3 or 4 or 5 or 51400.newchem");
     CommandArgument_Bool_OrDefault(PRINT_PATHS, False);
     CommandArgument_Bool_OrDefault(MAIN_LOGGING, False);
     CommandArgument_Bool_OrDefault(DETAILS, False);
     EndCommandArguments;

     // Process multiple phases if requested.

     if ( !PHASE.IsInt( ) )
     {    vec<int> phase;
          ParseIntSet( "{" + PHASE + "}", phase );
          for ( int i = 0; i < phase.isize( ); i++ )
          {    SystemSucceed( "Snorgle PHASE=" + ToString(phase[i])
                    + " INSTANCE=" + INSTANCE );    }
          Scram(0);    }

     int rev = 52330;

     if ( TRACE_RID >= 0 ) TRACE_PID = TRACE_RID/2;

     // String SAMPLE = "rhody";
     // String INSTANCE = "bugsy/rhody";

     String SAMPLE = "NA12878";

     String work_dir = "/wga/scr4/jaffe/GapToy/" + INSTANCE;

     // Phase 0.

     if ( PHASE == "0" )
     {  
          HyperBasevector hb;
          BinaryReader::readFile( work_dir + "/a.200/a.hbv", &hb );

          // Compute distance to end for each vertex.

          const int max_dist = 10000000; // dangerous!
          vec<int> D;
          cout << Date( ) << ": computing distances to end" << endl;
          int hbk = hb.K( );
          DistancesToEndFunc( hb, [hbk]( bvec const& bv ){ return bv.size()-hbk+1; },
                max_dist, True, D );
          BinaryWriter::writeFile( work_dir + "/dist_to_end", D );
          cout << Date( ) << ": done" << endl;
          Scram(0);    }

     // Phase 1.

     if ( PHASE == "1" )
     {    cout << Date( ) << ": loading" << endl;
          ReadPathVec paths( work_dir + "/a.200/a.paths" );
          int64_t npids = paths.size( ) / 2;
          vecbasevector bases( work_dir + "/data/frag_reads_orig.fastb" );
          int64_t N = bases.size( );
          VecPQVec pquals( work_dir + "/data/frag_reads_orig.qualp" );
          HyperBasevectorX hb;
          BinaryReader::readFile( work_dir + "/a.200/a.hbx", &hb );
          vec<int> inv;
          BinaryReader::readFile( work_dir + "/a.200/a.inv", &inv );
     
          // Heuristics.

          const int K = 40;
          const int64_t max_mult = 80;
          const int adaptn = 20;
          const int batches = 32;
          const int64_t bins = 20;
          const int qtop = 30;
          const int passes = 4;
     
          // Get distance to end for each vertex.

          vec<int> D;
          BinaryReader::readFile( work_dir + "/dist_to_end", &D );

          // Kill reads having adapter.
     
          double clock = WallClockTime( );
          cout << Date( ) << ": killing reads having adapter" << endl;
          String adapt 
               = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTATCTCGTATGCCGTCTTCTGCTTG";
          String adapt20 = adapt;
          adapt20.resize(adaptn);

          // Identify incompletely placed pairs, having some quality.  Ignore reads 
          // containing adapter.

          vec<Bool> un( bases.size( ), False );
          const int min_qual = 1000;
          const int max_prox = 500;
          #pragma omp parallel for
          for ( int64_t pid = 0; pid < npids; pid++ )
          {    int64_t id1 = 2*pid, id2 = 2*pid+1;
               const ReadPath &p1 = paths[id1], &p2 = paths[id2];

               if ( pid == TRACE_PID ) 
               {    cout << "checking PID = " << TRACE_PID << endl;
                    cout << "p1: " << p1.getOffset( )
                         << " : " << printSeq(p1) << endl;
                    cout << "p2: " << p2.getOffset( )
                         << " : " << printSeq(p2) << endl;    }

               // Test proximity.
     
               int d1 = -1, d2 = -1;
               if ( p1.size( ) > 0 )
               {    int start1 = p1.getOffset( );
                    int stop1 = start1 + bases[id1].isize( );
                    for ( int j = 1; j < (int) p1.size( ); j++ )
                         stop1 -= hb.Kmers( p1[j] );
                    int d1a = hb.Bases( p1.back( ) ) - stop1
                         + D[ hb.ToRight( p1.back( ) ) ];
                    int d1b = start1 + D[ hb.ToRight( inv[ p1.front( ) ] ) ];
                    d1 = Min( d1a, d1b );    }
               if ( p2.size( ) > 0 )
               {    int start2 = p2.getOffset( );
                    int stop2 = start2 + bases[id2].isize( );
                    for ( int j = 1; j < (int) p2.size( ); j++ )
                         stop2 -= hb.Kmers( p2[j] );
                    int d2a = hb.Bases( p2.back( ) ) - stop2
                         + D[ hb.ToRight( p2.back( ) ) ];
                    int d2b = start2 + D[ hb.ToRight( inv[ p2.front( ) ] ) ];
                    d2 = Min( d2a, d2b );    }
               if ( d1 > max_prox && d2 > max_prox ) continue;
               if ( d1 > max_prox && d2 == -1 ) continue;
               if ( d2 > max_prox && d1 == -1 ) continue;

               if ( pid == TRACE_PID ) cout << "past first check" << endl;

               /*
               if ( p1.size( ) > 0 && p2.size( ) > 0 )
               {    int e1 = p1.back( ), e2 = p2.back( );
                    // if ( hb.From( to_right[e1] ).nonempty( ) ) continue;
                    // if ( hb.From( to_right[e2] ).nonempty( ) ) continue;
                    if ( inv[e1] == e2 ) continue;    }
               if ( pid == TRACE_PID ) cout << "past second check" << endl;
               */

               int qsum1 = 0, qsum2 = 0;
               qvec qv1, qv2;
               pquals[id1].unpack( &qv1 );
               pquals[id2].unpack( &qv2 );
               for ( int j = 0; j < (int) qv1.size( ); j++ )
                    if ( qv1[j] > 2 ) qsum1 += qv1[j];
               for ( int j = 0; j < (int) qv2.size( ); j++ )
                    if ( qv2[j] > 2 ) qsum2 += qv2[j];

               if ( p1.size( ) == 0 && qsum1 < min_qual ) continue;
               if ( p2.size( ) == 0 && qsum2 < min_qual ) continue;
               if ( pid == TRACE_PID ) cout << "past third check" << endl;

               Bool have_ad = False;
               String s1 = bases[id1].ToString( );
               for ( int j = 0; j <= bases[id1].isize( ) - adaptn; j++ )
               {    if ( s1.Contains( adapt20, j ) )
                    {    have_ad = True;
                         break;    }    }
               if (have_ad) continue;
               String s2 = bases[id2].ToString( );
               for ( int j = 0; j <= bases[id2].isize( ) - adaptn; j++ )
               {    if ( s2.Contains( adapt20, j ) )
                    {    have_ad = True;
                         break;    }    }
               if (have_ad) continue;
               if ( pid == TRACE_PID ) cout << "past fourth check" << endl;
               un[id1] = un[id2] = True;    }

          // Build trash.

          vecbasevector trashb;
          trashb.reserve( 2 * Sum(un) );
          cout << Date( ) << ": building trash" << endl;
          int count = 0;
          vec<int> ids;
          for ( int64_t pid = 0; pid < npids; pid++ )
          {    int64_t id1 = 2*pid, id2 = 2*pid+1;
               if ( un[id1] )
               {    ids.push_back( id1, id2 );
                    trashb.push_back( bases[id1] );
                    trashb.push_back( bases[id2] );
                    count++;    }    }
          PRINT2( count, npids );
          BinaryWriter::writeFile( work_dir + "/hard_ids", ids );
          cout << PERCENT_RATIO( 3, count, npids ) << endl;
     
          // Reduce reads.
     
          vecqualvector qualsx( ids.size( ) );
          for ( int j = 0; j < ids.isize( ); j++ )
               pquals[ids[j]].unpack( &qualsx[j] );
          Destroy(pquals);
          vecbasevector basesx( ids.size( ) );
          for ( int j = 0; j < ids.isize( ); j++ )
               basesx[j] = bases[ ids[j] ];
          Destroy(bases);

          // Make kmers.

          cout << Date( ) << ": making kmers" << endl;
          cout << "current mem usage = " << MemUsageGBString( ) << endl;
          vec< triple<kmer<K>,int,int> > kmers_plus;
          MakeKmerLookup2x( trashb, kmers_plus );
          cout << "current mem usage = " << MemUsageGBString( ) << endl;
     
          // Set up to find friends.
     
          cout << Date( ) << ": finding friends" << endl;
          vec<int64_t> bstart(batches+1);
          for ( int64_t i = 0; i <= batches; i++ )
               bstart[i] = ( (int64_t) kmers_plus.size( ) * i ) / batches;
          #pragma omp parallel for schedule(dynamic, 1)
          for ( int64_t i = 1; i < batches; i++ )
          {    int64_t& s = bstart[i];
               while( s > 0 && kmers_plus[s].first == kmers_plus[s-1].first )
               {    s--;    }    }

          // Find friends.

          vec< quad<int64_t,int64_t,int,Bool> > xfriends;
          ReportPeakMem( );
          vec< vec< quad<int64_t,int64_t,int,Bool> > > friendsi(batches);
          #pragma omp parallel for schedule(dynamic, 1)
          for ( int64_t bi = 0; bi < batches; bi++ )
          {    vec< triple<kmer<K>,int,int> > left, right;
               for ( int64_t i = bstart[bi]; i < bstart[bi+1]; i++ )
               {    int64_t j;
                    for ( j = i + 1; j < bstart[bi+1]; j++ )
                         if ( kmers_plus[j].first != kmers_plus[i].first ) break;
                    if ( j - i == 1 ) continue;

                    // Handle the case where the stack is too large.

                    if ( j - i > max_mult )
                    {    
                         // Have to explore both right and left.  With a little
                         // reengineering, it would enough to create one copy.

                         right.resize( j - i );
                         left.resize( j - i );
                         memcpy( &right[0], &kmers_plus[i],
                              sizeof( triple<kmer<K>,int,int> ) * (j-i) );
                         memcpy( &left[0], &kmers_plus[i],
                              sizeof( triple<kmer<K>,int,int> ) * (j-i) );
                         ProcessRight( right, basesx, qualsx, max_mult, 
                              friendsi, bi, ids, TRACE_PID );
                         ProcessLeft( left, basesx, qualsx, max_mult,
                              friendsi, bi, ids, TRACE_PID ); 
                         i = j - 1;
                         continue;    }

                    // Process stack.

                    for ( int64_t k1 = i; k1 < j; k1++ )
                    for ( int64_t k2 = i; k2 < j; k2++ )
                    {    if ( k1 == k2 ) continue;
                         int id1 = kmers_plus[k1].second; 
                         int id2 = kmers_plus[k2].second;
                         if ( id1 == id2 ) continue;
                         int pos1 = kmers_plus[k1].third; 
                         int pos2 = kmers_plus[k2].third;
                         const basevector &b1 = trashb[id1], &b2 = trashb[id2];

                         // Reads may not agree on next left base.  This means
                         // left on the first read after making the first read fw
                         // if it is not already.

                         if ( pos1 >= 0 && pos2 >= 0 )
                         {    if ( pos1 > 0 && pos2 > 0 && b1[pos1-1] == b2[pos2-1] )
                                   continue;    }
                         else if ( pos1 < 0 && pos2 < 0 )
                         {    if ( -pos1-1 > 0 && -pos2-1 > 0
                                   && b1[-pos1-2] == b2[-pos2-2] ) 
                              {    continue;    }    }
                         else if ( pos1 >= 0 && pos2 < 0 )
                         {    if ( pos1 > 0 && -pos2+K <= b2.isize( )
                                   && b1[pos1-1] == 3 - b2[-pos2+K-1] )
                              {    continue;    }    }
                         else if ( pos1 < 0 && pos2 >= 0 )
                         {    if ( -pos1-1 > 0 && pos2 + K < b2.isize( )
                                   && 3 - b1[-pos1-2] == b2[pos2+K] )
                              {    continue;    }    }

                         // Trace.

                         if ( TRACE_PID >= 0 && TRACE_PID == ids[id1] / 2 )
                         {
                              #pragma omp critical
                              {    cout << "see overlap of " << ids[id1] << " with "
                                        << ids[id2] << endl;    }    }
                         if ( TRACE_PID >= 0 && TRACE_PID == ids[id2] / 2 )
                         {
                              #pragma omp critical
                              {    cout << "see overlap of " << ids[id2] << " with "
                                        << ids[id1] << endl;    }    }

                         // Form alignment.

                         Bool fw2 = !( pos1 < 0 ^ pos2 < 0 );
                         if ( pos1 < 0 ) pos1 = -pos1-1;
                         if ( pos2 < 0 ) pos2 = -pos2-1;
                         if ( !fw2 ) pos2 = trashb[id2].isize( ) - pos2 - K;
                         int offset = pos1 - pos2;

                         // Kill if hq diff, then save.

                         if ( HQDiff( id1, id2, offset, fw2, basesx, qualsx ) ) 
                              continue;
                         friendsi[bi].push( id1, id2, offset, fw2 );    }
                    i = j - 1;    }    
               UniqueSort( friendsi[bi] );    }

          // Collate friends.
     
          cout << Date( ) << ": pushing" << endl;
          ReportPeakMem( );
          vec<int64_t> starts( batches + 1, 0 );
          for ( int64_t bi = 0; bi < batches; bi++ )
               starts[bi+1] = starts[bi] + friendsi[bi].jsize( );
          xfriends.resize( starts.back( ) );
          #pragma omp parallel for
          for ( int64_t bi = 0; bi < batches; bi++ )
          {    if ( friendsi[bi].nonempty( ) )
               {    memcpy( &xfriends[ starts[bi] ], &friendsi[bi][0], 
                         sizeof( quad<int64_t,int64_t,int,Bool> ) 
                         * friendsi[bi].jsize( ) );    }
               Destroy( friendsi[bi] );    }

          int64_t x1 = xfriends.size( ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          ParallelUniqueSort( xfriends );
          int64_t x2 = xfriends.size( ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          cout << "reduced to: " << PERCENT_RATIO( 3, x2, x1 ) << endl; // XXXXXXXXX

          Destroy(kmers_plus), Destroy(trashb);
          if ( TRACE_PID >= 0 ) Scram(0);

          // Fix ids.

          cout << Date( ) << ": fixing ids" << endl;
          #pragma omp parallel for
          for ( int64_t i = 0; i < xfriends.jsize( ); i++ )
          {    xfriends[i].first = ids[ xfriends[i].first ];
               xfriends[i].second = ids[ xfriends[i].second ];    }



          /*
          Bool DUMP = True;
          if (DUMP)
          {
          for ( int64_t i = 0; i < xfriends.jsize( ); i++ ) // XXXXXXXXXXXXXXXXXXXXX
          {    int id1 = xfriends[i].first, id2 = xfriends[i].second; // XXXXXXXXXXX
               int offset = xfriends[i].third; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               int fw2 = xfriends[i].fourth; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               PRINT4( id1, id2, offset, fw2 );    

               // PATCHING BUG:
               if ( !fw2 )
               {    offset = offset;
                    // offset = (offset + rl2 - rl1);
                    swap( id1, id2 );
                    PRINT4( id1, id2, offset, fw2 );    }

                    } // XXXXXXXXXXXXXXXXXXXXXXXXXXXX
          }
          */



/*
             pos1                        **************
         ================================>
             <=========================================
             pos2

offset = - (offset + rl2 - rl1)

offset = pos1 - pos2
*/


          cout << "\n";
          int64_t nfriends = xfriends.size( );
          PRINT(nfriends);
          cout << Date( ) << ": " << TimeSince(clock) << " used" << endl;
          ReportPeakMem( );
          cout << endl;
     
          // Save friends.
     
          cout << Date( ) << ": saving friends" << endl;
          BinaryWriter::writeFile( work_dir + "/xfriends", xfriends );

          cout << Date( ) << ": creating index" << endl;
          vec<int64_t> index( N+1, -1 );
          for ( int64_t i = xfriends.jsize( ) - 1; i >= 0; i-- )
               index[ xfriends[i].first ] = i;
          index[N] = xfriends.size( );
          for ( int64_t i = N-1; i >= 0; i-- )
               if ( index[i] < 0 ) index[i] = index[i+1];
          cout << Date( ) << ": writing index" << endl;
          BinaryWriter::writeFile( work_dir + "/xfriends.index", index );

          ReportPeakMem( );
          cout << Date( ) << ": done" << endl;
          Scram(0);    }

     // Phase 2.

     if ( PHASE == "2" )
     {    double clock = WallClockTime( );

          Bool FP_LOGGING = False;
          Bool PRINT_FRIENDS = False;
          Bool PRINT_STACKS = False;
          Bool PRINT_JSTACK = False;

          if ( TRACE_PID >= 0 )
          {    MAIN_LOGGING = True;
               FP_LOGGING = True;
               PRINT_FRIENDS = True;
               PRINT_STACKS = True;
               PRINT_JSTACK = True;    }

          vec< quad<int64_t,int64_t,int,Bool> > xfriends;
          cout << Date( ) << ": reading friends" << endl;
          BinaryReader::readFile( work_dir + "/xfriends", &xfriends );

          /*
          VirtualMasterVec<basevector> basesx_main
               ( work_dir + "/data/frag_reads_orig.fastb" );
          */

          MasterVec<basevector> basesx( work_dir + "/data/frag_reads_orig.fastb" );

          // VirtualMasterVec<qualvector> qualsx_main
          //      ( work_dir + "/data/frag_reads_orig.qualb" );

          // Need to virtualize access:
          VecPQVec pqualsx( work_dir + "/data/frag_reads_orig.qualp" );

          // Set up bins and closures.

          const int64_t bins = 400;
          vec<int64_t> starts( bins + 1 );
          for ( int j = 0; j <= bins; j++ )
               starts[j] = ( j * xfriends.jsize( ) ) / bins;
          for ( int j = bins - 1; j > 0; j-- )
          {    while( starts[j] > 0 && xfriends[ starts[j] ].first/2 
                    == xfriends[ starts[j]-1 ].first/2 )
               {    starts[j]--;    }    }
          vec<vec<basevector>> closures(bins);
          vec<vec<int64_t>> closure_ids(bins);

          // Parallel loop over each of the bins m.

          cout << Date( ) << ": starting main loop" << endl;
          #pragma omp parallel for schedule(dynamic,1)
          for ( int64_t m = 0; m < bins; m++ )
          {    int64_t start = starts[m], stop = starts[m+1];
               // VirtualMasterVec<basevector> basesx(basesx_main);
               int64_t last_id = -1;
               readstack last_stack, stack, jstack;
               PairsManager pairs;
               vecbasevector bases;
               vecqualvector quals;
               Friends aligns;
               vec<Bool> suspect;
               basevector last_con, con, jcon;
               vec<basevector> p;
               vec<int> ids, off;
               for ( int64_t i = start; i < stop; i++ )
               {    // Reminder: xfriends[i] is data (id1,id2,offset,fw2).
                    int id = xfriends[i].first;
                    if ( TRACE_PID >= 0 && id/2 != TRACE_PID ) continue;

                    // Restrict to those reads for which id1 = id.

                    int64_t j;
                    for ( j = i + 1; j < stop; j++ )
                         if ( xfriends[j].first != id ) break;

                    if (MAIN_LOGGING)
                    {    cout << "\n" << xfriends[i].first << " has " << j - i 
                              << " friends" << endl;    }
          
                    // if ( j - i <= 20 )
                    if (PRINT_FRIENDS)
                    {    for ( int64_t k = i; k < j; k++ )
                         {    cout << ( xfriends[k].fourth ? "+" : "-" )
                                   << xfriends[k].second 
                                   << "." << xfriends[k].third << endl;    }    }
     
                    // Load reads.

                    ids = {id};
                    for ( int64_t k = i; k < j; k++ )
                         ids.push_back( xfriends[k].second );
                    double lclock = WallClockTime( );
                    bases.clear( );
                    quals.clear( );

                    for ( int k = 0; k < ids.isize( ); k++ )
                    {    bases.push_back( basesx[ ids[k] ] );
                         qvec q;
                         pqualsx[ids[k]].unpack( &q );
                         quals.push_back(q);    }
                    if (MAIN_LOGGING)
                         cout << TimeSince(lclock) << " used loading reads" << endl;

                    // Cap quality scores.

                    vec<Bool> done( quals.size( ), False );
                    CapQualityScores( quals, done );

                    // Build read stack.

                    double rclock = WallClockTime( );
                    aligns.clear( );
                    for ( int64_t k = i; k < j; k++ )
                    {    aligns.push_back( Friend( k - i + 1, xfriends[k].third, 
                              !xfriends[k].fourth ) );    }
                    stack.Initialize( 0, aligns, 0, aligns.size( ), 
                         readstack::right_extended, bases, quals, pairs, False );
                    for ( int64_t k = i; k < j; k++ )
                    {    stack.SetId( k-i, ids[k-i] );
                         stack.SetPid( k-i, ids[k-i]/2 );    }
                    stack.Raise1(0);
                    stack.MotifDiff(1,suspect);
                    stack.Erase(suspect);

                    const int q_solid = 30;
                    stack.HighQualDiff( q_solid, 1, suspect );
                    if (MAIN_LOGGING) PRINT2( int(suspect[0]), int(suspect[1]) );
                    stack.Erase(suspect);
     
                    if ( xfriends[i].first % 2 == 1 ) stack.Reverse( );
     
                    if (PRINT_STACKS) stack.Print(cout);
     
                    // Look for overlaps.
          
                    con = stack.Consensus1( );
                    if ( xfriends[i].first % 2 == 1 
                         && xfriends[i].first == last_id + 1 )
                    {    double oclock = WallClockTime( );
                         off = GetOffsets1( last_stack, stack, 0, 0 );
                         if (MAIN_LOGGING)
                         {    cout << TimeSince(oclock) << " used getting offsets" 
                                   << endl;
                              cout << "off = " << printSeq(off) << endl;    }
                         const int max_offsets = 2;
                         if ( TRACE_PID >= 0 
                              && ( off.empty( ) || off.isize( ) > max_offsets ) )
                         {    cout << "\n";
                              con.Print( cout, "consensus" );
                              basevector conl = last_stack.Consensus1( );
                              conl.Print( cout, "last_consensus" );    
                              cout << "\n";    }

                         if ( off.isize( ) <= max_offsets )
                         {    for ( int o = 0; o < off.isize( ); o++ )
                              {    if (MAIN_LOGGING)
                                   {    cout << "\nusing offset " << off[o] 
                                             << endl;    }
                                   jstack = last_stack;
                                   jstack.Merge( stack, off[o] );

                                   // What does this do and does it assume
                                   // reads have same length?

                                   if ( -2*off[o] >= jstack.Cols( ) ) continue;

                                   // Trim so that stack does not extend beyond
                                   // the originating fragment.

                                   int left = 0, right = jstack.Cols( );
                                   if ( off[o] < 0 ) left = -off[o];
                                   int ext_right = last_stack.Cols( ) -
                                        ( off[o] + stack.Cols( ) );
                                   if ( ext_right > 0 ) right -= ext_right;
                                   jstack.Trim( left, right );

                                   // Ignore very short fragments.

                                   if ( jstack.Cols( ) < 100 ) continue;

                                   jstack.SortByPid( 
                                        jstack.Pid(0), 0, last_stack.Rows( ) );
                                   jstack.Unique( );
                                   // jstack.Raise1(0), jstack.Raise1(1);
                                   jstack.HighQualDiff( q_solid, 2, suspect );
                                   // if ( suspect[0] || suspect[1] ) continue;
                                   if (MAIN_LOGGING)
                                   {    PRINT2( int(suspect[0]), int(suspect[1]) );
                                        PRINT( jstack.Rows( ) );    }
                                   jstack.Erase(suspect);
                                   if (MAIN_LOGGING) PRINT( jstack.Rows( ) );
               
                                   // Clean stack.
          
                                   for ( int j = 0; j < jstack.Cols( ); j++ )
                                   {    vec<int> qsum( 4, 0 ); 
                                        vec<int> ids( 4, vec<int>::IDENTITY );
                                        for ( int l = 0; l < jstack.Rows( ); l++ )
                                        {    if ( jstack.Def( l, j ) )
                                             {    qsum[ jstack.Base( l, j ) ]
                                                       += jstack.Qual( l, j );    
                                                       }    }
                                        ReverseSortSync( qsum, ids );
                                        if ( qsum[0] >= 100 
                                             && qsum[0] >= 10 * qsum[1]
                                             && qsum[1] < 100 )
                                        {    for ( int l = 0; l < jstack.Rows( ); 
                                                  l++ )
                                             {    if ( jstack.Def( l, j ) 
                                                       && jstack.Base( l, j ) 
                                                            != ids[0] )
                                                  {    jstack.SetBase(l, j, ids[0]);
                                                       jstack.SetQual( l, j, 0 );    
                                                            }    }    }    }
     
                                   // Print stack.
     
                                   if (PRINT_JSTACK) 
                                   {    cout << "\njoint stack:\n";
                                        jstack.Print(cout);    }
                                   jcon = jstack.Consensus1( );
                                   if (MAIN_LOGGING) 
                                        jcon.Print( cout, "joint_consensus" );    
                                   FindPaths( jstack, p, FP_LOGGING );
                                   const int max_paths = 5;
                                   if ( p.isize( ) <= max_paths )
                                   {    for ( int c = 0; c < p.isize( ); c++ )
                                        {    closures[m].push_back( p[c] );    
                                             closure_ids[m].push_back(id);    }    }
                                   if ( MAIN_LOGGING || PRINT_PATHS )
                                   {    
                                        #pragma omp critical
                                        {    cout << "id = " << id << ", off[0] = "
                                                  << off[0] << ", found " 
                                                  << p.size( ) << " paths";
                                             if ( p.isize( ) > max_paths )
                                                  cout << " (too many)";
                                             cout << "\n";
                                             if ( p.size( ) <= 12 )
                                             {    for ( int j = 0; j < p.isize( ); 
                                                       j++ )
                                                  {    p[j].Print( cout, 
                                                            "p" + ToString(j+1) );
                                                       }    }    }    }    }    }

                         /*
                         basevector con1 = last_con;
                         basevector con2 = con;
                         vec<int> offsets;
                         const int L = 20;
                         vecbasevector cons;
                         cons.push_back(con1);
                         cons.push_back(con2);
                         vec< triple<kmer<L>,int,int> > kmers_plus;
                         MakeKmerLookup0( cons, kmers_plus );
                         for ( int m = 0; m < kmers_plus.isize( ); m++ )
                         {    int k;
                              for ( k = m + 1; k < kmers_plus.isize( ); k++ )
                              {    if ( kmers_plus[k].first != kmers_plus[m].first ) 
                                        break;    }
                              for ( int l1 = m; l1 < k; l1++ )
                              for ( int l2 = m; l2 < k; l2++ )
                              {    if ( kmers_plus[l1].second != 0 ) continue;
                                   if ( kmers_plus[l2].second != 1 ) continue;
                                   int o = kmers_plus[l1].third 
                                        - kmers_plus[l2].third;
                                   offsets.push_back(o);    }
                              m = k - 1;    }
                         UniqueSort(offsets);
                         cout << "\noffsets = " << printSeq(offsets) << "\n" << endl;
     
                         // Score the offsets.
     
                         for ( int m = 0; m < offsets.isize( ); m++ )
                         {    vec<int> perfs;
                              int count = 0;
                              for ( int p2 = 0; p2 < con2.isize( ); p2++ )
                              {    int p1 = p2 + offsets[m];
                                   if ( p1 < 0 ) continue;
                                   if ( p1 >= con1.isize( ) ) break;
                                   if ( con1[p1] != con2[p2] )
                                   {    if ( count >= L ) perfs.push_back(count);
                                        count = 0;    }
                                   else count++;    }
                              if ( count >= L ) perfs.push_back(count);
                              ReverseSort(perfs);
                              if (MAIN_LOGGING)
                              {    cout << "offset = " << offsets[m] << ", perfs = "
                                        << printSeq(perfs) << endl;    }    }
                         if (MAIN_LOGGING) cout << "\n";     
                         */
     
                         }

                    if (MAIN_LOGGING) cout << TimeSince(rclock) << " used" << endl;
     
                    last_id = xfriends[i].first;
                    last_con = con;
                    last_stack = stack;
                    i = j - 1;    }    }

          if ( TRACE_PID >= 0 ) Scram(0);
          cout << Date( ) << ": merging closures" << endl;
          vecbasevector closuresx;
          vec<int64_t> closure_idsx;
          for ( int m = 0; m < closures.isize( ); m++ )
          for ( int j = 0; j < closures[m].isize( ); j++ )
          {    closuresx.push_back( closures[m][j] );
               closure_idsx.push_back( closure_ids[m][j] );    }
          closuresx.WriteAll( work_dir + "/closures.fastb" );
          BinaryWriter::writeFile( work_dir + "/closures.ids", closure_idsx );
          cout << Date( ) << ": found " << closuresx.size( ) << " closures" << endl;

          // Double closures.
     
          vecbasevector clo2(closuresx);
          clo2.Append(closuresx);
          for ( int64_t i = 0; i < (int64_t) closuresx.size( ); i++ )
               clo2[ closuresx.size( ) + i ].ReverseComplement( );
          vec<int64_t> closure_idsx2(closure_idsx);
          closure_idsx2.append(closure_idsx);

          // Build K=100 lookup table for genome.

          const int K = 100;
          vecbasevector genome( work_dir + "/a.200/a.fastb" );
          vecbasevector gx(genome);
          HyperBasevectorX hb;
          BinaryReader::readFile( work_dir + "/a.200/a.hbx", &hb );
          vec<int> cstart( hb.E( ) ), cstop( hb.E( ) );
          for ( int e = 0; e < hb.E( ); e++ )
          {    int Kbig = hb.K( );
               const basevector& u = hb.EdgeObject(e);
               cstart[e] = 0, cstop[e] = u.isize( ) - Kbig + K;
               int v = hb.ToLeft(e), w = hb.ToRight(e);
               if ( hb.To(v).size( ) == 1 )
               {    if ( Kbig - K < cstop[e] ) cstart[e] = Kbig - K;    }
               if ( hb.To(w).size( ) == 1 ) cstop[e] = u.isize( );
               gx[e].SetToSubOf( gx[e], cstart[e], cstop[e] - cstart[e] );    }

          // Find closures having a unique genomic start/stop.  And kill closures 
          // that exactly match the assembly.

          cout << Date( ) << ": finding closures having unique start/stop" << endl;
          vec< pair<int,int> > astart( clo2.size( ), make_pair(-1,0) );
          vec< pair<int,int> > astop( clo2.size( ), make_pair(-1,0) );

          // Kill closures having a unique K-mer (K smaller).

          {    const int K = 80;
               vecbasevector gc(genome);
               gc.Append(clo2);
               vec< triple<kmer<K>,int,int> > kmers_plus;
               MakeKmerLookup0( gc, kmers_plus );
               for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
               {    int64_t j;
                    for ( j = i + 1; j < kmers_plus.jsize( ); j++ )
                         if ( kmers_plus[j].first != kmers_plus[i].first ) break;
                    if ( j - i == 1 )
                    {    int id = kmers_plus[i].second;
                         if ( id >= (int) genome.size( ) )
                              clo2[ id - genome.size( ) ].resize(0);    }
                    i = j - 1;    }    }

          // Build genome lookup.

          cout << Date( ) << ": building genome lookup" << endl;
          vec< triple<kmer<K>,int,int> > gkmers_plus;
          MakeKmerLookup0( gx, gkmers_plus );

          // Find starts and stops.

          cout << Date( ) << ": finding start/stop" << endl;
          #pragma omp parallel for
          for ( int i = 0; i < (int) clo2.size( ); i++ )
          {    basevector& c = clo2[i];
               if ( c.isize( ) < K ) continue;
               kmer<K> x1, x2;
               x1.SetToSubOf( c, 0 );
               x2.SetToSubOf( c, c.isize( ) - K );
               int64_t low1 = LowerBound1( gkmers_plus, x1 );
               int64_t high1 = UpperBound1( gkmers_plus, x1 );
               if ( high1 - low1 == 1 )
               {    int gid = gkmers_plus[low1].second;
                    int gstart = gkmers_plus[low1].third + cstart[gid];

                    // Kill matching closures.

                    // if ( gstart + c.isize( ) <= genome[gid].isize( ) )
                    {    int e = gid, p = gstart + K - 1, pos1;
                         for ( pos1 = K; pos1 < c.isize( ); pos1++ )
                         {    ++p;
                              if ( p < hb.Bases(e) ) 
                              {    if ( hb.O(e)[p] != c[pos1] ) break;    }
                              else
                              {    Bool match = False;
                                   int v = hb.ToRight(e);
                                   for ( int r = 0; 
                                        r < (int) hb.From(v).size( ); r++ )
                                   {    int en = hb.IFrom(v,r);
                                        if ( hb.O(en)[hb.K( )-1] == c[pos1] )
                                        {    e = en;
                                             p = hb.K( )-1;
                                             match = True;
                                             break;    }    }
                                   if ( !match ) break;    }    }
                         if ( pos1 == c.isize( ) )
                         {    c.resize(0);
                              continue;    }    }

                    // Mark start.

                    astart[i] = make_pair( gid, gstart );
                    if (DETAILS)
                    {    
                         #pragma omp critical
                         {    cout << "pair " << closure_idsx2[i]/2 << ", closure " 
                                   << i << " starts at " << gid
                                   << "." << gstart << endl;    }    }    }
               int64_t low2 = LowerBound1( gkmers_plus, x2 );
               int64_t high2 = UpperBound1( gkmers_plus, x2 );

               // Test for consistent stops.  Note asymmetry: should do this for
               // starts too.

               const int max_stops = 5;
               if ( high2 - low2 > 1 && high2 - low2 <= max_stops )
               {    Bool consistent = True;
                    for ( int64_t i1 = low2; i1 < high2; i1++ )
                    for ( int64_t i2 = i1+1; i2 < high2; i2++ )
                    {    int e1 = gkmers_plus[i1].second;
                         int e2 = gkmers_plus[i2].second;
                         int offset = gkmers_plus[i1].third - gkmers_plus[i2].third;
                         for ( int p1 = 0; p1 < genome[e1].isize( ); p1++ )
                         {    int p2 = p1 - offset;
                              if ( p2 < 0 || p2 >= genome[e2].isize( ) ) continue;
                              if ( genome[e1][p1] != genome[e2][p2] )
                              {    consistent = False;
                                   goto next;    }    }    }
                    next:
                    if (consistent)
                    {    vec<int> right( high2 - low2 );
                         right[0] = genome[ gkmers_plus[low2].second ].size( );
                         for ( int64_t i = low2 + 1; i < high2; i++ )
                         {    int offset 
                                   = gkmers_plus[i].third - gkmers_plus[low2].third;
                              right[i-low2] = offset 
                                   + genome[ gkmers_plus[i].second ].size( );    }
                         int R = Max(right);
                         int r;
                         for ( r = 0; r < right.isize( ); r++ )
                              if ( right[r] == R ) break;
                         low2 += r;
                         high2 = low2 + 1;    }    }

               // Now handle unique placements.

               if ( high2 - low2 == 1 )
               {    int gid = gkmers_plus[low2].second;
                    int gstart = gkmers_plus[low2].third + cstart[gid];
                    astop[i] = make_pair( gid, gstart + K );
                    if (DETAILS)
                    {    
                         #pragma omp critical
                         {    cout << "pair " << closure_idsx2[i]/2 << ", closure " 
                                   << i << " stops at " << gid
                                   << "." << gstart + K << endl;    }    }    }    }
          BinaryWriter::writeFile( work_dir + "/closures.astart", astart );
          BinaryWriter::writeFile( work_dir + "/closures.astop", astop );
          clo2.WriteAll( work_dir + "/closures2.fastb" );

          // Build lookup table for closures.
     
          cout << Date( ) << ": making lookup" << endl;
          vec< triple<kmer<K>,int,int> > kmers_plus;
          MakeKmerLookup0( clo2, kmers_plus );
     
          // Find overlaps between closures.
     
          cout << Date( ) << ": looking for overlaps" << endl;
          vec< vec< pair<int,int> > > O( clo2.size( ) );
          int64_t overs = 0;
          const int max_overs = 200;
          for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
          {    int64_t j;
               for ( j = i + 1; j < kmers_plus.jsize( ); j++ )
                    if ( kmers_plus[j].first != kmers_plus[i].first ) break;
               if ( j - i > max_overs )
               {    i = j - 1;
                    continue;    }
               for ( int64_t k1 = i; k1 < j; k1++ )
               {    if ( kmers_plus[k1].third != 0 ) continue;
                    for ( int64_t k2 = i; k2 < j; k2++ )
                    {    if ( k2 == k1 ) continue;
                         int id1 = kmers_plus[k1].second; 
                         int id2 = kmers_plus[k2].second;
                         Bool mismatch = False;
                         for ( int p1 = K; p1 < clo2[id1].isize( ); p1++ )
                         {    int p2 = kmers_plus[k2].third + p1;
                              if ( p2 >= clo2[id2].isize( ) ) break;
                              if ( clo2[id1][p1] != clo2[id2][p2] )
                              {    mismatch = True;
                                   break;    }    }
                         if (mismatch) continue;
                         O[id2].push( id1, kmers_plus[k2].third );
                         overs++;    }    }
               i = j - 1;    }
          cout << Date( ) << ": there are " << ToStringAddCommas(overs)
               << " overlaps" << endl;
          BinaryWriter::writeFile( work_dir + "/closures.overlaps", O );

          // Done.

          cout << Date( ) << ": done" << endl;
          cout << "\n" << TimeSince(clock) << " used" << endl;
          ReportPeakMem( );
          Scram(0);    }

     // Phase 3.

     if ( PHASE == "3" )
     {    cout << Date( ) << ": loading" << endl;
          vecbasevector closures( work_dir + "/closures.fastb" );

          HyperBasevector hb;
          BinaryReader::readFile( work_dir + "/a.200/a.hbv", &hb );
          vec<int> to_right;
          hb.ToRight(to_right);

          vecbasevector edges;
          for ( int e = 0; e < hb.E( ); e++ )
               edges.push_back( hb.EdgeObject(e) );
          const int K = 60;
          const int K2 = hb.K( );
          vec< triple<kmer<K>,int,int> > kmers_plus;
     
          cout << Date( ) << ": making kmer lookup" << endl;
          MakeKmerLookup0( edges, kmers_plus );
          vecbasevector patches;
          cout << "\n" << Date( ) << ": looking up closures" << endl;
          for ( int i = 0; i < (int) closures.size( ); i++ )
          {    const basevector& C = closures[i];
               if ( C.isize( ) >= K )
               {    vec<int> p1, p2;

                    // Look for the first and last kmers of C.

                    kmer<K> c1, c2;
                    c1.SetToSubOf( C, 0 );
                    c2.SetToSubOf( C, C.isize( ) - K );
                    int64_t low1 = LowerBound1( kmers_plus, c1 );
                    int64_t high1 = UpperBound1( kmers_plus, c1 );
                    int64_t low2 = LowerBound1( kmers_plus, c2 );
                    int64_t high2 = UpperBound1( kmers_plus, c2 );
                    for ( int64_t j = low1; j < high1; j++ )
                         p1.push_back( kmers_plus[j].second );
                    for ( int64_t j = low2; j < high2; j++ )
                         p2.push_back( kmers_plus[j].second );
     
                    // Ignore perfect placements.
     
                    Bool perf = False;
                    for ( int64_t j = low1; j < high1; j++ )
                    {    int e = kmers_plus[j].second, epos = kmers_plus[j].third;
                         int en = hb.EdgeLengthBases(e);
                         int cpos = 0, cn = C.size( );
                         ForceAssertEq( C[cpos], hb.EdgeObject(e)[epos] ); // XXXXXX
                         while ( cpos < cn - K )
                         {    if ( epos < en - K )
                              {    if ( hb.EdgeObject(e)[epos+K] == C[cpos+K] )
                                   {    epos++;
                                        cpos++;    }
                                   else break;    }
                              else // epos = en - K
                              {    int v = to_right[e];
                                   Bool found = False;
                                   for ( int m = 0; m < hb.From(v).isize( ); m++ )
                                   {    int f = hb.IFrom( v, m );
                                        {    if ( hb.EdgeObject(f)[K2-1] 
                                                  == C[cpos+K] )
                                             {    e = f;
                                                  en = hb.EdgeLengthBases(e);
                                                  epos = K2 - K;
                                                  cpos++;
                                                  found = True;
                                                  break;    }    }    }
                                   if ( !found ) break;    }    }
                         if ( cpos == cn - K )
                         {    perf = True;
                              break;    }    }    
                    if (perf) continue;

                    // Print.

                    const int max_patches = 10;
                    if ( p1.isize( ) <= max_patches && p2.isize( ) <= max_patches
                         && p1.isize( ) * p2.isize( ) <= max_patches && p1 != p2 )
                    {    if (MAIN_LOGGING)
                         {    cout << "\n" << i << ": " << printSeq(p1) << " --> " 
                                   << printSeq(p2) << endl;    }    }
                    else continue;
                    if (MAIN_LOGGING) C.Print( cout, "closure" );

                    // Glue.

                    for ( int64_t m1 = low1; m1 < high1; m1++ )
                    for ( int64_t m2 = low2; m2 < high2; m2++ )
                    {    int e1 = kmers_plus[m1].second;
                         int ne1 = hb.EdgeLengthBases(e1);
                         int epos1 = kmers_plus[m1].third;
                         int e2 = kmers_plus[m2].second;
                         int ne2 = hb.EdgeLengthBases(e2);
                         int epos2 = kmers_plus[m2].third;
                         int K2p = K2 + 1;
                         if ( epos1 + K >= K2p && ne2 - epos2 >= K2p )
                         {    const basevector& E1 = hb.EdgeObject(e1);
                              const basevector& E2 = hb.EdgeObject(e2);
                              basevector R( E2, epos2, K2p );
                              basevector patch( E1, epos1 + K - K2p, K2p );
                              if ( C.isize( ) >= 2*K )
                              {    patch = Cat( patch, basevector( 
                                        C, K, C.isize( ) - 2*K ), R );    }
                              else
                              {    int over = 2*K - C.isize( );
                                   patch.resize( patch.isize( ) - over );
                                   patch = Cat( patch, R );    }
                              if (MAIN_LOGGING) patch.Print( cout, "patch" );    
                              patches.push_back(patch);    }    }    }    }
     
          cout << "\n" << Date( ) << ": writing " << patches.size( )
               << " patches" << endl;
          patches.WriteAll( work_dir + "/new_stuff.exp" );
          cout << "\n" << Date( ) << ": done" << endl;
          Scram(0);    }

     // Phase 4.

     if ( PHASE == "4" )
     {    String I = INSTANCE;
          if ( I.Contains( "/" ) ) I = I.Before( "/" );
          SystemSucceed( "GapToy." + ToString(rev) + " X=all SAMPLE=" + SAMPLE +
               " INSTANCE=" + INSTANCE 
               + " START_WITH_PATCHES=True SAVE_PATCHED=True EVALUATE=False BETSYBOB=True FIN=3" 
               );    }    }
