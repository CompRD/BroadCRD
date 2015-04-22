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
#include "FetchReads.h"
#include "PrintAlignment.h"
#include "Qualvector.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "pairwise_aligners/ClusterAligner.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "pairwise_aligners/SmithWatBandedQual.h"
#include "paths/LongReadTools.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "paths/long/ValidateRef.h"
#include "polymorphism/Edit.h"

void ProcessAlignment( const int pass, const int id, const look_align& la, 
     const basevector& b, const qualvector& q, const basevector& r, int& qual_sum, 
     int& max_perf, const int verbosity, const int trace_id, ostream& outx )
{    
     const align& az = la.a;
     vec<ho_interval> perfs2;
     az.PerfectIntervals2( b, r, perfs2 );
     for ( int j = 0; j < perfs2.isize( ); j++ )
          max_perf = Max( max_perf, perfs2[j].Length( ) );
     int p1 = az.pos1( ), p2 = az.pos2( );
     for ( int j = 0; j < az.Nblocks( ); j++ ) 
     {    if ( az.Gaps(j) > 0 ) // deletion
          {    int qx = q[p1];
               if ( p1 > 0 && qx > q[p1-1] ) qx = q[p1-1];
               qual_sum += az.Gaps(j) * qx;
               p2 += az.Gaps(j);    }
          if ( az.Gaps(j) < 0 ) // insertion
          {    for ( int z = 0; z < -az.Gaps(j); z++ )
                    qual_sum += q[p1+z];
               p1 -= az.Gaps(j);    }
          for ( int x = 0; x < az.Lengths(j); x++ ) 
          {    if ( b[p1] != r[p2] ) qual_sum += q[p1];
               ++p1; ++p2;    }    }
     if ( verbosity >= 3 )
     {    
          #pragma omp critical
          {    outx << "read " << id << " aligns to " << ( pass == 1 ? "A" : "R" ) 
                    << " with " << la.Errors( ) << " errors" << " and qual sum = " 
                    << qual_sum << ", start = " << az.pos2( ) << endl;     }    }
     if ( id == trace_id || verbosity >= 4 )
     {    
          #pragma omp critical
          {    outx << "\nalignment of read " << id << endl;
               PrintVisualAlignment( True, outx, b, r, az, q );    }    }    }

void ProcessRead( const int id, const vecbasevector& R, const vecbasevector& A, 
     const int K, const VecIntPairVec& Rlocs,
     const VecIntPairVec& Alocs, const vecbasevector& bases,
     const vecqualvector& quals, const int min_perfect_match, const int min_qual,
     vec<int>& terrs, const int verbosity, const int trace_id, ostream& out )
{    
     terrs.resize_and_set( 2, 0 );
     vec<int> errs( 2, 1000000000 );
     int max_perf = 0;
     ostringstream outx;
     for ( int pass = 1; pass <= 2; pass++ )
     {    const vecbasevector& P = ( pass == 1 ? A : R );
          const VecIntPairVec& locs = ( pass == 1 ? Alocs : Rlocs );
          vec<look_align> aligns;
          ClusterAligner( bases[id], P, K, locs, aligns );
          for ( int h = 0; h < aligns.isize( ); h++ )
          {    basevector b = bases[id];
               qualvector q = quals[id];
               if ( aligns[h].Rc1( ) )
               {    b.ReverseComplement( );
                    q.ReverseMe( );    }

               // Realign, taking account of quality scores.

               vec<int> offsets;
               align& az = aligns[h].a;
               int p1 = az.pos1( ), p2 = az.pos2( );
               offsets.push_back( p1 - p2 );
               for ( int j = 0; j < az.Nblocks( ); j++ ) 
               {    if ( az.Gaps(j) > 0 )  
                    {    p2 += az.Gaps(j);
                         offsets.push_back( p1 - p2 );    }
                    if ( az.Gaps(j) < 0 ) 
                    {    p1 -= az.Gaps(j);
                         offsets.push_back( p1 - p2 );    }
                    p1 += az.Lengths(j); p2 += az.Lengths(j);    }
               int mino = Min(offsets), maxo = Max(offsets);
               int offset = (mino + maxo)/2;
               int bandwidth = Max( offset - mino, maxo - offset ) + 2;
               int errors;
               SmithWatBandedQual( b, P[0], q, offset, bandwidth, 
                    aligns[h].a, errors, 0, 1, 1 );

               // Flatten low quality scores.

               for ( int j = 0; j < (int) q.size( ); j++ )
                    if ( q[j] < min_qual ) q[j] = 0;

               // Process alignment.

               int qual_sum = 0;
               ProcessAlignment( pass, id, aligns[h], b, q, P[0], qual_sum,
                    max_perf, verbosity, trace_id, outx );
               errs[pass-1] = Min( errs[pass-1], qual_sum );    }    }
     const int max_min_score = 100;
     if ( errs[0] != errs[1] && max_perf >= min_perfect_match
          && Min( errs[0], errs[1] ) <= max_min_score )
     {    
          #pragma omp critical
          {    out << outx.str( );
               if ( verbosity >= 2 ) PRINT3_TO( out, id, errs[0], errs[1] );    
               if ( errs[0] < errs[1] ) terrs[1] = errs[1] - errs[0];
               else terrs[0] = errs[0] - errs[1];    }    }    }

// CompareRefs: compare two 'references' R and A by aligning reads (bases, quals)
// to them.  As implemented, assumes only one record in each reference.

void CompareRefs( const vecbasevector& R, const vecbasevector& A, const int K, 
     const vecbasevector& bases, const vecqualvector& quals, const int min_support,
     const int min_perfect_match, const int min_qual, vec<int>& total_errs, 
     const int verbosity, const int trace_id, ostream& out )
{
     // Create indices to facilitate alignment.

     VecIntPairVec Rlocs, Alocs;
     CreateGlocs( R, K, Rlocs );
     CreateGlocs( A, K, Alocs );

     // Align the reads to the edited reference.
     
     vec<int> incr;
     incr.push_back( 256, 128, 64, 32, 16, 8, 4, 2, 1 );
     total_errs.resize_and_set( 2, 0 );
     for ( int pass = 0; pass < incr.isize( ); pass++ )
     {    
          #pragma omp parallel for
          for ( int id = 0; id < (int) bases.size( ); id += incr[pass] )
          {    Bool done = False;
               for ( int j = 0; j < pass; j++ )
                    if ( id % incr[j] == 0 ) done = True;
               if (done) continue;
               vec<int> terrs;
               ProcessRead( id, R, A, K, Rlocs, Alocs, bases, quals, 
                    min_perfect_match, min_qual, terrs, verbosity, trace_id, out );
               #pragma omp critical
               {    for ( int j = 0; j < 2; j++ )
                         total_errs[j] += terrs[j];    }    }
          if ( total_errs[1] >= min_support ) break;    }    }

void DefineDiffs( const basevector& a, const basevector& r,
     const int mismatch, const int gap_open, const int gap_extend,
     vec< pair<int,edit0> >& edits, const int verbosity )
{
     alignment al;
     SmithWatAffineParallel( a, r, al, true, true, mismatch, gap_open, gap_extend );
     if ( verbosity >= 1 )
     {    cout << "\nalignment of assembly to reference:\n";
          PrintVisualAlignment( True, cout, a, r, al );    }
     align alx(al);
     edits.clear( );
     int p1 = alx.pos1( ), p2 = alx.pos2( );
     for ( int j = 0; j < alx.Nblocks( ); j++ ) 
     {    if ( alx.Gaps(j) > 0 )  
          {    edits.push( p2, edit0( DELETION, alx.Gaps(j) ) );
               p2 += alx.Gaps(j);    }
          if ( alx.Gaps(j) < 0 ) 
          {    basevector b( a, p1, -alx.Gaps(j) );
               edits.push( p2, edit0( INSERTION, b.ToString( ) ) );
               p1 -= alx.Gaps(j);    }
          for ( int x = 0; x < alx.Lengths(j); x++ ) 
          {    if ( a[p1] != r[p2] )
                    edits.push( p2, edit0( SUBSTITUTION, as_base( a[p1] ) ) );
               ++p1; ++p2;    }    }    }

void Compare( const basevector& a, const basevector& r, const vecbasevector& bases, 
     const vecqualvector& quals, const int K, const int min_ratio, 
     const int min_support, const int min_perfect_match, const int max_gap,
     const int min_qual, const int mismatch, const int gap_open, 
     const int gap_extend, const int verbosity, const int VALIDATE_TRACE_ID,
     const Bool VALIDATE_DISPLAY_ALL )
{
     // Align the assembly to the reference.

     vec< pair<int,edit0> > edits;
     DefineDiffs( a, r, mismatch, gap_open, gap_extend, edits, verbosity );

     // Go through the edits.

     cout << "testing " << edits.size( ) << " edits" << endl;
     Bool succeed = True;
     for ( int ei = 0; ei < edits.isize( ); ei++ )
     {    
          // Find nearby edits.

          int ej;
          for ( ej = ei + 1; ej < edits.isize( ); ej++ )
               if ( edits[ej].first - edits[ej-1].first > max_gap ) break;

          // Set up references.

          vecbasevector R(1);
          R[0] = r;
          vecbasevector A(R);

          // Go through edits.

          for ( int ek = ej - 1; ek >= ei; ek-- )
          {
               // Define edit.
     
               const edit0& e = edits[ek].second;
               int pos = edits[ek].first;

               // Introduce edits.

               if ( e.etype == DELETION )
               {    basevector x;
                    for ( int j = 0; j < pos; j++ )
                         x.push_back( A[0][j] );
                    for ( int j = pos + e.n; j < A[0].isize( ); j++ )
                         x.push_back( A[0][j] );
                    A[0] = x;    }
               if ( e.etype == INSERTION )
               {    basevector x;
                    for ( int j = 0; j < pos; j++ )
                         x.push_back( A[0][j] );
                    for ( int j = 0; j < e.seq.isize( ); j++ )
                         x.push_back( as_char( e.seq[j] ) );
                    for ( int j = pos; j < A[0].isize( ); j++ )
                         x.push_back( A[0][j] );
                    A[0] = x;    }
               if ( e.etype == SUBSTITUTION )
                    A[0].Set( pos, as_char( e.seq[0] ) );    }

          // Compare references.

          vec<int> total_errs;
          ostringstream out;
          CompareRefs( R, A, K, bases, quals, min_support, min_perfect_match,
               min_qual, total_errs, verbosity, VALIDATE_TRACE_ID, out );
          Bool OK = ( total_errs[1] >= min_support 
               && total_errs[1] >= min_ratio * total_errs[0] );
          if ( !OK ) succeed = False;
          if ( !OK || VALIDATE_DISPLAY_ALL )
          {    cout << "\n";
               for ( int ek = ei; ek < ej; ek++ )
               {    const edit0& e = edits[ek].second;
                    int pos = edits[ek].first;
                    cout << "testing edit " << ek << ": " << e << " at " << pos 
                         << "\n";    }
               cout << out.str( );
               cout << "alt errors = " << total_errs[0]
                    << ", ref errors = " << total_errs[1] << endl;    }
          ei = ej - 1;    }
     if (succeed) cout << "\nvalidation succeeded\n";
     else cout << "\nvalidation failed\n";    }

void ValidateRef( const SupportedHyperBasevector& shb, const String& TMP,
     const long_logging_control& log_control, const long_logging& logc )
{    
     // Set up.

     Bool verbosity = logc.verb[ "VALIDATE" ];
     double clock = WallClockTime( );
     cout << "\n";

     // Define heuristics.

     const int min_ratio = 5;           // alt qual sum must be this times lower
     const int min_support = 250;       // use enough reads until see this qual sum
     const int min_perfect_match = 40;  // minimum perfect match for read alignments
     const int mismatch = 2;            // mismatch penalty for ref alignment
     const int gap_open = 3;            // gap open penalty for ref alignment
     const int gap_extend = 1;          // gap extend penalty for ref alignment
     const int K = 11;                  // seed size for read alignments
     const int max_gap = 20;            // max separation to cluster edits
     const int min_qual = 3;            // zero out quality scores below this

     // Load data.

     int edge_id = logc.VALIDATE;
     ForceAssert( edge_id < shb.EdgeObjectCount( ) );
     vecbasevector A;
     A.push_back( shb.EdgeObject(edge_id) );
     vecbasevector RR( *(log_control.G) );
     vecbasevector bases( TMP + "/0:0.fastb" );
     bases.ReadAll( TMP + "/1:0.fastb", True );
     vecqualvector quals( TMP + "/0:0.qualb" );
     quals.ReadAll( TMP + "/1:0.qualb", True );

     // Try a couple tag sizes.

     Bool found_match = False;
     for ( int tag = 100; tag >= 80; tag -= 20 )
     {    if (found_match) break;

          // Sanity check data.

          if( A[0].isize( ) < 2 * tag )
          {    cout << "Assembly edge too short to validate.\n";
               cout << "Validation failed.\n\n";
               return;    }

          // Go through reference contigs.

          for ( int r = 0; r < (int) RR.size( ); r++ )
          {    if (found_match) break;
               vecbasevector R;
               R.push_back( RR[r] );
               cout << "comparing to reference contig " << r << endl;
               if( R[0].isize( ) < 2 * tag ) 
               {    cout << "it is too short\n";
                    continue;    }

               // Two passes because assembly might be in either of two orientations.

               for ( int pass = 1; pass <= 2; pass++ )
               {    if ( pass == 1 ) cout << "attempting forward validation\n";
                    else cout << "\nattempting reverse validation\n";
                    if ( pass ==  2 ) A[0].ReverseComplement( );
                    basevector r1( R[0], 0, tag ); 
                    basevector r2( R[0], R[0].isize( ) - tag, tag );
                    basevector a1( A[0], 0, tag ); 
                    basevector a2( A[0], A[0].isize( ) - tag, tag );
                    String As = A[0].ToString( );
                    String a1s = a1.ToString( ), a2s = a2.ToString( );
                    String Rs = R[0].ToString( );
                    String r1s = r1.ToString( ), r2s = r2.ToString( );
                    vec<int> a1locs, a2locs, r1locs, r2locs;
                    for ( int p = 0; p <= A[0].isize( ); p++ )
                    {    if ( As.Contains( r1s, p ) ) 
                         {    cout << "found r1 at " << p << endl;
                              r1locs.push_back(p);    }
                         if ( As.Contains( r2s, p ) ) 
                         {    cout << "found r2 at " << p 
                                   << " of " << A[0].isize( ) << endl;
                              r2locs.push_back(p);    }    }
                    for ( int p = 0; p <= R[0].isize( ); p++ )
                    {    if ( Rs.Contains( a1s, p ) ) 
                         {    cout << "found a1 at " << p << endl;
                              a1locs.push_back(p);    }
                         if ( Rs.Contains( a2s, p ) ) 
                         {    cout << "found a2 at " << p 
                                   << " of " << R[0].isize( ) << endl;
                              a2locs.push_back(p);    }    }    
     
                    // Attempt to find common ends for assembly and reference.
     
                    Bool match = False;
                    basevector a, r;
                    if ( r1locs.solo( ) && r2locs.solo( ) && r1locs[0] < r2locs[0] )
                    {    a = basevector( 
                              A[0], r1locs[0], r2locs[0] + tag - r1locs[0] );
                         r = R[0];
                         match = True;    }
                    else if ( a1locs.solo( ) && a2locs.solo( ) 
                         && a1locs[0] < a2locs[0] )
                    {    r = basevector( 
                              R[0], a1locs[0], a2locs[0] + tag - a1locs[0] );
                         a = A[0];
                         match = True;    }
                    else if ( a1locs.solo( ) && r2locs.solo( ) )
                    {    r = basevector( 
                              R[0], a1locs[0], R[0].isize( ) - a1locs[0] );
                         a = basevector( A[0], 0, r2locs[0] + tag );
                         match = True;    }
                    else if ( r1locs.solo( ) && a2locs.solo( ) )
                    {    r = basevector( R[0], 0, a2locs[0] + tag );
                         a = basevector( 
                              A[0], r1locs[0], A[0].isize( ) - r1locs[0] );
                         match = True;    }
                    if (match)
                    {    found_match = True;
                         Compare( a, r, bases, quals, K, min_ratio, min_support, 
                              min_perfect_match, max_gap, min_qual, mismatch, 
                              gap_open, gap_extend, verbosity, 
                              logc.VALIDATE_TRACE_ID, logc.VALIDATE_DISPLAY_ALL );
                         break;    }    }    }    }

     cout << "\n";
     REPORT_TIME( clock, "used validating" );    }
