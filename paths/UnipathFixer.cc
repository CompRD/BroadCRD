///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// UnipathFixer.  Attempt to fix unipaths by aligning the reads back to them and
// using these alignments to edit putative defects in the unipaths, for example
// by making note of poorly covered loci.  The alignments place the uncorrected
// reads, with the goal of finding the best home(s) for each, even though those
// homes have errors.  The method for finding the best home is heuristic and almost
// certainly wrong sometimes, but perhaps good enough.  The algorithm generates
// a stack of bases and quality scores for each position on a unipath, which are
// then used as input to the editing steps.
//
// UNDER CONSTRUCTION!

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "FeudalMimic.h"
#include "Intvector.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "Qualvector.h"
#include "VecUtilities.h"
#include "feudal/BinaryStream.h"
#include "graph/Digraph.h"
#include "math/Functions.h"
#include "math/Matrix.h"
#include "paths/GetNexts.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerBasket.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/UnibaseUtils.h"
#include "paths/Unipath.h"
#include "paths/UnipathFixerTools.h"
#include "util/SearchFastb2Core.h"
#include "PairsManager.h"

class int5 {
     public:
     unsigned int x[5];
     friend Bool operator==( const int5& I1, const int5& I2 )
     {    return memcmp( I1.x, I2.x, 20 ) == 0;    }
     friend Bool operator!=( const int5& I1, const int5& I2 )
     {    return memcmp( I1.x, I2.x, 20 ) != 0;    }
     friend Bool operator<( const int5& I1, const int5& I2 )
     {    return memcmp( I1.x, I2.x, 20 ) < 0;    }
};

// Compare a triple by its second entry, then its third, then its first.

Bool cmp_middle( const triple<int64_t,int64_t,int>& T1, 
     const triple<int64_t,int64_t,int>& T2 )
{    if ( T1.second < T2.second ) return True;
     if ( T1.second > T2.second ) return False;
     if ( T1.third < T2.third ) return True;
     if ( T1.third > T2.third ) return False;
     if ( T1.first < T2.first ) return True;
     return False;    }

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_Int(K);
     CommandArgument_String_OrDefault(IN_HEAD, "all_reads");
     CommandArgument_String_OrDefault(OUT_HEAD, "all_reads.fixed");
     CommandArgument_Bool_OrDefault(PRINT_ALIGNMENTS, False);
     CommandArgument_Bool_OrDefault(PRINT_SEGMENTS, False);
     CommandArgument_Bool_OrDefault(SHOW_ALL, False);
     CommandArgument_Int_OrDefault(MAX_PLACEMENTS, 50);
     CommandArgument_Bool_OrDefault(VALIDATE, False);
     CommandArgument_Bool_OrDefault(WRITE, True);
     CommandArgument_Bool_OrDefault(FISHING_VERBOSE, False);
     CommandArgument_Bool_OrDefault(DUMP_MISSING_GENOMIC_KMERS, False);
     CommandArgument_Bool_OrDefault(PRINT_ASSEMBLIES, False);
     CommandArgument_String_OrDefault(PRINT_ASSEMBLIES_VERBOSE_LIST, "");
     CommandArgument_Bool_OrDefault(USE_JUMPS, True);
     EndCommandArguments;

     // Start.

     cout << Date( ) << ": start" << endl;
     double clock = WallClockTime( );
     String data_dir = PRE + "/" + DATA, run_dir = PRE + "/" + DATA + "/" + RUN;
     Mkdir777( run_dir + "/tmp" );
     String temp_file = run_dir + "/tmp/UnipathFixer.fastb";
     String temp_file2 = run_dir + "/tmp/UnipathFixer.2.fastb";
     String unifile = run_dir + "/tmp/UnipathFixer.unibases.fastb";
     vecbasevector genome;
     if (VALIDATE) genome.ReadAll( data_dir + "/genome.fastb" );
     if (VALIDATE) ForceAssertEq( K, 80 );
     vec<int> print_assemblies_verbose_list;
     ParseIntSet( PRINT_ASSEMBLIES_VERBOSE_LIST, print_assemblies_verbose_list );

     // Load unibases.

     cout << Date( ) << ": load unibases" << endl;
     vecbasevector unibases( run_dir + "/" + IN_HEAD + ".unibases.k80" );

     // Compute summary stats from initial unibases.

     vec<basevector> old_missing_genomic_kmers;
     int64_t old_non_genomic_kmers = 0, old_N50_unibase = 0;
     if (VALIDATE)
     {    cout << Date( ) << ": computing summary stats from old unibases" << endl;
          if ( K == 80 )
          {    UnibaseSummaryStats<80>( unibases, genome, old_missing_genomic_kmers,
                    old_non_genomic_kmers, old_N50_unibase );    }
          else ForceAssert( 0 == 1 );    }

     // Load unibases and remove hanging ends.

     cout << Date( ) << ": removing hanging ends" << endl;
     {    vecKmerPath paths, pathsrc, unipaths;
          vec<tagged_rpint> pathsdb, unipathsdb;
          ReadsToPathsCoreY( unibases, K, paths, pathsrc, pathsdb );
          Unipath( paths, pathsrc, pathsdb, unipaths, unipathsdb );
          digraph A;
          BuildUnipathAdjacencyGraph( paths, pathsrc, pathsdb, unipaths,
               unipathsdb, A );
          KmerBaseBroker( K, paths, pathsrc, pathsdb, unibases );
          HyperKmerPath h;
          BuildUnipathAdjacencyHyperKmerPath( K, A, unipaths, h );
          RemoveHangingEnds( h, &KmerPath::KmerCount, 250, 5.0 );
          KmerBaseBroker kbb( K, paths, pathsrc, pathsdb, unibases );
          unibases.clear( ); 
          for ( int i = 0; i < h.EdgeObjectCount( ); i++ )
               unibases.push_back( kbb.Seq( h.EdgeObject(i) ) );
          unibases.WriteAll(unifile);    }
     cout << "\nnew unibases:\n";
     for ( size_t i = 0; i < unibases.size( ); i++ )
          unibases[i].Print( cout, "unibases_" + ToString(i) );
     cout << "\n";

     // Set up ancillary data structures for unibases.

     cout << Date( ) << ": defining ancillary data structures for unibases" << endl;
     size_t nuni = unibases.size( );
     vec< vec<int> > nexts;
     GetNexts( K, unibases, nexts );
     cout << Date( ) << ": nexts computed" << endl;
     vec<int> to_rc;
     UnibaseInvolution( unibases, to_rc );

     // Fish for connections between dead ends.

     cout << Date( ) << ": looking for dead ends" << endl;
     vecbasevector dead_ends;
     vec<int> dead_fw, dead_rc;
     for ( size_t u = 0; u < nuni; u++ )
     {    if ( nexts[u].empty( ) ) dead_fw.push_back(u);
          if ( nexts[ to_rc[u] ].empty( ) ) dead_rc.push_back(u);    }
     const int KS = 40;
     const int KS_MAX_PLACEMENTS = 100;
     const int KS_K = 40;
     for ( int i = 0; i < dead_fw.isize( ); i++ )
     {    int u = dead_fw[i];
          basevector e;
          e.SetToSubOf( unibases[u], unibases[u].isize( ) - KS, KS );
          dead_ends.push_back(e);    }
     for ( int i = 0; i < dead_rc.isize( ); i++ )
     {    int u = dead_rc[i];
          basevector e;
          e.SetToSubOf( unibases[u], 0, KS );
          dead_ends.push_back(e);    }
     dead_ends.WriteAll(temp_file);
     vec<vecbasevector> extenders( dead_ends.size( ) );
     vec<vecqualvector> extenders_q( dead_ends.size( ) );
     int passes = ( USE_JUMPS ? 2 : 1 );
     for ( int pass = 1; pass <= passes; pass++ )
     {    String fastb = run_dir + "/" + ( pass == 1 ? "frag" : "jump" ) 
               + "_reads_filt.fastb";
          String qualb = run_dir + "/" + ( pass == 1 ? "frag" : "jump" ) 
               + "_reads_filt.qualb";
          vec< triple<int64_t,int64_t,int> > xaligns;
          vecbasevector reads(fastb);
          vecqualvector quals(qualb);
          SearchFastb2( temp_file, fastb, KS_K, &xaligns, 0, KS_MAX_PLACEMENTS );
          for ( int j = 0; j < xaligns.isize( ); j++ )
          {    if ( xaligns[j].third >= 0 )
               {    extenders[ xaligns[j].first ].push_back( 
                         reads[ xaligns[j].second ] );
                    extenders_q[ xaligns[j].first ].push_back( 
                         quals[ xaligns[j].second ] );    }
               else
               {    basevector b;
                    b.ReverseComplement( reads[ xaligns[j].second ] );
                    extenders[ xaligns[j].first ].push_back(b);
                    extenders_q[ xaligns[j].first ].push_back(
                         Reverse( quals[ xaligns[j].second ] ) );    }    }    }
     cout << Date( ) << ": found extenders" << endl;
     vecbasevector hedges;
     for ( int i = 0; i < dead_fw.isize( ); i++ )
     {    int u = dead_fw[i];
          if (FISHING_VERBOSE)
          {    cout << "\nforward from u = " << u << endl;
               unibases[u].Print( cout, "u=" + ToString(u) );
               for ( size_t j = 0; j < extenders[i].size( ); j++ )
                    extenders[i][j].Print( cout, j );    }
          HyperBasevector h;
          KmerBasket( extenders[i], extenders_q[i], h, 40, 10, 0, cout );
          for ( int j = 0; j < h.EdgeObjectCount( ); j++ )
               hedges.push_back( h.EdgeObject(j) );
          if (FISHING_VERBOSE) cout << h;    }
     for ( int i = 0; i < dead_rc.isize( ); i++ )
     {    int u = dead_rc[i], ip = i + dead_fw.isize( );
          if (FISHING_VERBOSE)
          {    cout << "\nbackward from u = " << u << endl;
               unibases[u].Print( cout, "u=" + ToString(u) );
               for ( size_t j = 0; j < extenders[ip].size( ); j++ )
                    extenders[ip][j].Print( cout, j );    }
          HyperBasevector h;
          KmerBasket( extenders[ip], extenders_q[ip], h, 40, 10, 0, cout );
          for ( int j = 0; j < h.EdgeObjectCount( ); j++ )
               hedges.push_back( h.EdgeObject(j) );
          if (FISHING_VERBOSE) cout << h;    }
     hedges.WriteAll(temp_file);
     vec< triple<int64_t,int64_t,int> > yaligns;
     SearchFastb2( temp_file, run_dir + "/" + IN_HEAD + ".unibases.k80", 
          K, &yaligns );
     vec<Bool> to_deletex( hedges.size( ), False );
     for ( int i = 0; i < yaligns.isize( ); i++ )
          to_deletex[ yaligns[i].first ] = True;
     for ( size_t i = 0; i < hedges.size( ); i++ )
          if ( hedges[i].isize( ) <= K ) to_deletex[i] = True;
     dead_ends.WriteAll(temp_file);
     vecbasevector new_stuff;
     for ( size_t i = 0; i < hedges.size( ); i++ )
          if ( !to_deletex[i] ) new_stuff.push_back( hedges[i] );
     new_stuff.WriteAll(temp_file2);
     SearchFastb2( temp_file, temp_file2, KS_K, &yaligns );
     vec< vec<int> > new_from( new_stuff.size( ) ), new_to( new_stuff.size( ) );
     for ( int j = 0; j < yaligns.isize( ); j++ )
     {    if ( yaligns[j].third >= 0 )
          {    if ( yaligns[j].first < dead_fw.isize( ) )
               {    new_from[ yaligns[j].second ].push_back( 
                         dead_fw[ yaligns[j].first ] );    }
               else 
               {    new_to[ yaligns[j].second ].push_back( 
                         dead_rc[ yaligns[j].first 
                              - dead_fw.isize( ) ] );    }    }    }
     vecbasevector good_new_stuff;
     for ( size_t i = 0; i < new_stuff.size( ); i++ )
     {    if ( !new_from[i].solo( ) || !new_to[i].solo( ) ) continue;
          if ( unibases[ new_from[i][0] ].isize( ) < 2*K ) continue;
          if ( unibases[ new_to[i][0] ].isize( ) < 2*K ) continue;
          good_new_stuff.push_back( new_stuff[i] );    }

     /*
     cout << "\nnew stuff:\n";
     for ( size_t i = 0; i < new_stuff.size( ); i++ )
          new_stuff[i].Print( cout, i );
     for ( int j = 0; j < yaligns.isize( ); j++ )
     {    if ( yaligns[j].third >= 0 )
          {    cout << yaligns[j].second << " reaches ";
               if ( yaligns[j].first < dead_fw.isize( ) )
                    cout << "from " << dead_fw[ yaligns[j].first ] << endl;
               else 
               {    cout << " to " << dead_rc[ yaligns[j].first - dead_fw.isize( ) ]
                         << endl;    }    }    }
     */

     // return 0; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     
     // Define *some* heuristic constants.

     const int min_to_keep_with_holes = 200;

     // Align reads. 
     
     vec< triple<int64_t,int64_t,int> > ALIGNS, JALIGNS;
     AlignReadsToUnipaths( run_dir, "jump_reads_filt", "frag_reads_filt",
          "frag_reads_edit", USE_JUMPS, MAX_PLACEMENTS, unifile, ALIGNS, JALIGNS );

     // Trim the jumps.  First sort the alignments by the unibase they start on.
     // This is denovo trimming: the idea is that if a jump read is supported by
     // another jump read (up to some point), then we don't trim it (before that
     // point).  What we actually do is look at the pile of jump reads over a given
     // unibase, and walk across it using a sliding window whose size varies as a
     // function of the depth of reads at that point. NOTE THAT THIS TENDS TO TRIM
     // BACK THE JUMP TO BEFORE ITS FIRST ERROR.

     vecbasevector jreads;
     size_t njreads = 0;
     vec<int> trim_to, starts;
     if (USE_JUMPS)
     {    jreads.ReadAll( run_dir + "/jump_reads_filt.fastb" );
          njreads = jreads.size( );
          trim_to.resize( njreads, 0 );
          cout << Date( ) << ": start jump trimming" << endl;
          sort( JALIGNS.begin( ), JALIGNS.end( ), cmp_middle );
          starts.push_back(0);
          for ( size_t i = 0; i < JALIGNS.size( ); i++ )
          {    size_t j;
               for ( j = i + 1; j < JALIGNS.size( ); j++ )
                    if ( JALIGNS[j].second != JALIGNS[i].second ) break;
               starts.push_back(j);    
               i = j - 1;    }
          #pragma omp parallel for
          for ( size_t aa = 0; aa < starts.size( ) - 1; aa++ )
          {    size_t i = starts[aa], j = starts[aa+1];
               int u = JALIGNS[i].second;
               for ( j = i + 1; j < JALIGNS.size( ); j++ )
                    if ( JALIGNS[j].second != u ) break;
               int n = j - i, nu = unibases[u].size( ), maxread = 0;
               for ( int l = 0; l < n; l++ )
                    maxread = Max( maxread, jreads[ JALIGNS[i+l].first ].isize( ) );

               // Divide the unibase into chunks.  This is for efficiency.  In 
               // principle it does not change the answer but in practice it does 
               // have a small effect.

               const int uchunk = 400;
               int read_start = 0;
               for ( int chunk_start = 0; chunk_start < nu; 
                    chunk_start += uchunk - maxread )
               {    int chunk_stop = Min( nu, chunk_start + uchunk );
                    int chunk_len = chunk_stop - chunk_start;
                    while( read_start < n && JALIGNS[i+read_start].third 
                         < chunk_start )
                    {    read_start++;    }
                    int read_stop = read_start;
                    while( read_stop < n && JALIGNS[i+read_stop].third < chunk_stop )
                         read_stop++;
                    int nr = read_stop - read_start;

                    // Create a matrix having nr rows (one per read) and 
                    // chunk_len + maxread (columns).  Stuff in the reads, filling 
                    // empty space with blanks.  Then carry out "de novo" trimming:
                    // we don't use the unipaths except to anchor the reads with 
                    // respect to each other.  Roughly, the first part of a read is
                    // 'validated' by another read if they agree.

                    vec< vec<char> > M(nr);
                    vec< vec<Bool> > validated(nr);
                    for ( int l = 0; l < nr; l++ )
                    {    M[l].resize( chunk_len + maxread, ' ' );
                         validated[l].resize( chunk_len + maxread, False );
                         int rid = JALIGNS[i+read_start+l].first; 
                         int pos = JALIGNS[i+read_start+l].third - chunk_start;
                         for ( int r = 0; r < jreads[rid].isize( ); r++ )
                              M[l][pos+r] = as_base( jreads[rid][r] );    }
                    for ( int z = 0; z < chunk_len + maxread; z++ )
                    {    
                         // Count the number of bases in column z.
     
                         int nbases = 0;
                         for ( int l = 0; l < nr; l++ )
                              if ( M[l][z] != ' ' ) nbases++;

                         // Now form the matrix sub consisting of the 
                         // ceil( log_4(nbases) ) + 2 columns starting with column z.

                         if ( nbases == 0 ) continue;
                         int log4 = 0, nb = nbases;
                         while ( nb > 0 )
                         {    log4++;
                              nb = nb >> 2;    }
                         int cols = log4 + 2;
                         if ( z + cols > chunk_len + maxread ) continue;
                         vec< vec<char> > sub(nr);
                         for ( int l = 0; l < nr; l++ )
                         {    sub[l].resize(cols);
                              for ( int r = 0; r < cols; r++ )
                                   sub[l][r] = M[l][z+r];    }

                         // Sort the rows of sub.  Then look for identical rows. 
                         // These validate entries in column z.  Note that duplicate
                         // pairs would foul this test.  Have these been removed
                         // already?

                         vec<int> ids( nr, vec<int>::IDENTITY );
                         SortSync( sub, ids );
                         for ( int l1 = 0; l1 < nr; l1++ )
                         {    int l2;
                              for ( l2 = l1 + 1; l2 < nr; l2++ )
                                   if ( sub[l2] != sub[l1] ) break;
                              Bool good = True;
                              for ( int s = 0; s < cols; s++ )
                                   if ( sub[l1][s] == ' ' ) good = False;
                              if ( good && l2 - l1 > 1 )
                              {    for ( int r = l1; r < l2; r++ )
                                        validated[ ids[r] ][z] = True;    }
                              l1 = l2 - 1;    }    }

                    // Update the read trim points.

                    for ( int l = 0; l < nr; l++ )
                    {    int rid = JALIGNS[i+read_start+l].first; 
                         int pos = JALIGNS[i+read_start+l].third - chunk_start;
                         for ( int r = 0; r < jreads[rid].isize( ); r++ )
                         {    if ( !validated[l][pos+r] )
                              {    
                                   #pragma omp critical
                                   {    trim_to[rid] = Max( trim_to[rid], r );    }
                                   break;    }    }    }

                    // Update read_start.

                    read_start = read_stop;
                    while( read_start >= 1 && JALIGNS[i+read_start-1].third 
                         >= chunk_start + uchunk - maxread )
                    {    read_start--;    }
                    read_start--;    }    }    }

     // Set up data structure to track calls at each unibase base.

     VecIntVec calls[4], qcalls[4];
     for ( int i = 0; i < 4; i++ )
     {    Mimic( unibases, calls[i] ), Mimic( unibases, qcalls[i] );    }

     // Extend alignments through the unipath graph.

     int64_t total = 0;  // total number of bases aligned to unibases
     int64_t qtotal = 0; // total quality aligned
     int64_t qgood = 0;  // total quality aligned correctly
     cout << Date( ) << ": extending fragments alignments" << endl;
     Sort(ALIGNS);
     vecbasevector reads( run_dir + "/frag_reads_filt.fastb" );
     vecqualvector quals( run_dir + "/frag_reads_filt.qualb" );
     vec<segalign> SEGS;
     ExtendAligns( K, ALIGNS, reads, quals, unibases, nexts, PRINT_ALIGNMENTS,
          PRINT_SEGMENTS, calls, qcalls, total, qtotal, qgood, SEGS, False );
     vec<segalign> JSEGS;
     vecqualvector jquals;
     if (USE_JUMPS)
     {    cout << Date( ) << ": extending jump alignments" << endl;
          Sort(JALIGNS);
          jquals.ReadAll( run_dir + "/jump_reads_filt.qualb" );
          for ( size_t i = 0; i < njreads; i++ )
          {    if ( trim_to[i] < 50 ) trim_to[i] = 0;
               jreads[i].resize( trim_to[i] ), jquals[i].resize( trim_to[i] );    }
          ExtendAligns( K, JALIGNS, jreads, jquals, unibases, nexts, 
               PRINT_ALIGNMENTS, PRINT_SEGMENTS, calls, qcalls, total, qtotal, 
               qgood, JSEGS, False );    }

     // If two unipaths differ at just one position, and the support for one at that
     // position is way better than the support for the other, kill the weak one.

     cout << Date( ) << ": look for unipaths having a much better partner" << endl;
     vec<Bool> to_delete(nuni, False);
     vecbasevector unibases_trunc(unibases);
     for ( size_t i = 0; i < nuni; i++ )
          unibases_trunc[i].resize(20);
     unibases_trunc.WriteAll(temp_file);
     vec< triple<int64_t,int64_t,int> > aligns;
     SearchFastb2( temp_file, unifile, 20, &aligns, 0, MAX_PLACEMENTS );
     vec<Bool> delete_close_to_better( nuni, False );
     for ( size_t i = 0; i < aligns.size( ); i++ )
     {    int u1 = aligns[i].first, u2 = aligns[i].second;
          if ( unibases[u1].size( ) != unibases[u2].size( ) ) continue;
          if ( aligns[i].third != 0 ) continue;
          vec<int> diffs;
          for ( int j = 0; j < unibases[u1].isize( ); j++ )
               if ( unibases[u1][j] != unibases[u2][j] ) diffs.push_back(j);
          if ( !diffs.solo( ) || delete_close_to_better[u2] ) continue;
          int j = diffs[0];
          int nu1 = unibases[u1].size( ), ru1 = to_rc[u1];
          int nu2 = unibases[u2].size( ), ru2 = to_rc[u2];
          int D1 = unibases[u1][j], D2 = unibases[u2][j];
          int Q1 = qcalls[D1][u1][j] + qcalls[3-D1][ru1][nu1-j-1];
          int Q2 = qcalls[D2][u2][j] + qcalls[3-D2][ru2][nu2-j-1];
          if ( Q1 < 20 && Q2 >= 5 * Q1 ) 
          {    delete_close_to_better[u1] = True;
               to_delete[u1] = True;    }    }

     // Look for discrepancies and form the natural read piles over them, then
     // look for new kmers. 
     // IN PROGRESS!!!

     cout << Date( ) << ": looking for clusters" << endl;
     vec<int> X(4), Q(4);

     Sort(SEGS);
     vec<size_t> SEGS_START(nuni+1);
     {    size_t SEG_POS = 0;
          for ( size_t u = 0; u <= nuni; u++ )
          {    while( SEG_POS < SEGS.size( ) && SEGS[SEG_POS].u < u ) ++SEG_POS;
               SEGS_START[u] = SEG_POS;    }    }
     vec<size_t> JSEGS_START(nuni+1);
     if (USE_JUMPS)
     {    Sort(JSEGS);
          {    size_t SEG_POS = 0;
               for ( size_t u = 0; u <= nuni; u++ )
               {    while( SEG_POS < JSEGS.size( ) && JSEGS[SEG_POS].u < u ) 
                         ++SEG_POS;
                    JSEGS_START[u] = SEG_POS;    }    }    }
     reads.ReadAll( run_dir + "/frag_reads_filt.fastb" );
     if (USE_JUMPS)
     {    jreads.ReadAll( run_dir + "/jump_reads_filt.fastb" );
          for ( size_t i = 0; i < njreads; i++ )
               jreads[i].resize( trim_to[i] );    }
     for ( size_t u = 0; u < nuni; u++ )
     {    int nu = unibases[u].size( ), ru = to_rc[u];
          for ( int j = 0; j < nu; j++ )
          {    for ( int k = 0; k < 4; k++ )
                    Q[k] = qcalls[k][u][j] + qcalls[3-k][ru][nu-j-1];
               int D = unibases[u][j];
               if ( Q[D] == Max(Q) ) continue;

               // Find the reads passing through u.j.

               vec<segalign> nhood;
               vec<String> ids;
               vec<Bool> is_frag;
               for ( size_t l = SEGS_START[u]; l < SEGS_START[u+1]; l++ )
               {    int rid = SEGS[l].rid;
                    if ( SEGS[l].upos + reads[rid].isize( ) - SEGS[l].rpos <= j ) 
                         continue;
                    if ( SEGS[l].upos >= j ) break;
                    nhood.push_back( SEGS[l] );
                    ids.push_back( ToString(rid) + ".frag.fw" );
                    is_frag.push_back(True);   }
               if (USE_JUMPS)
               {    for ( size_t l = JSEGS_START[u]; l < JSEGS_START[u+1]; l++ )
                    {    int rid = JSEGS[l].rid;
                         if ( JSEGS[l].upos + jreads[rid].isize( ) 
                              - JSEGS[l].rpos <= j ) 
                         {    continue;    }
                         if ( JSEGS[l].upos >= j ) break;
                         nhood.push_back( JSEGS[l] );
                         ids.push_back( ToString(rid) + ".jump.fw" );
                         is_frag.push_back(False);   }    }

               // Get the reads from the reverse complement.

               int jr = nu - j - 1;
               for ( size_t l = SEGS_START[ru]; l < SEGS_START[ru+1]; l++ )
               {    int rid = SEGS[l].rid;
                    if ( SEGS[l].upos + reads[rid].isize( ) <= jr ) continue;
                    if ( SEGS[l].upos >= jr ) break;
                    int rpos = SEGS[l].rpos, upos = SEGS[l].upos, rpos_new = 0;
                    int upos_new = nu - ( upos + reads[rid].size( ) - rpos );
                    while( upos_new < 0 )
                    {    ++rpos_new;
                         ++upos_new;    }
                    if ( upos_new > j ) continue;
                    nhood.push( False, rid, rpos_new, u, upos_new );
                    ids.push_back( ToString(rid) + ".frag.rc" );
                    is_frag.push_back(True);    }
               if (USE_JUMPS)
               {    for ( size_t l = JSEGS_START[ru]; l < JSEGS_START[ru+1]; l++ )
                    {    int rid = JSEGS[l].rid;
                         if ( JSEGS[l].upos + jreads[rid].isize( ) <= jr ) continue;
                         if ( JSEGS[l].upos >= jr ) break;
                         int rpos = JSEGS[l].rpos, upos = JSEGS[l].upos; 
                         int rpos_new = 0;
                         int upos_new = nu - ( upos + jreads[rid].size( ) - rpos );
                         while( upos_new < 0 )
                         {    ++rpos_new;
                              ++upos_new;    }
                         if ( upos_new > j ) continue;
                         nhood.push( False, rid, rpos_new, u, upos_new );
                         ids.push_back( ToString(rid) + ".jump.rc" );
                         is_frag.push_back(False);    }    }
               if ( nhood.empty( ) ) continue;

               // Display the matrix associated to these reads.

               int ustart = nhood.front( ).upos, ustop = 0;
               for ( int l = 0; l < nhood.isize( ); l++ )
               {    int rid = nhood[l].rid;
                    const basevector& R = ( is_frag[l] ? reads[rid] : jreads[rid] );
                    ustop = Max( ustop, 
                         nhood[l].upos + R.isize( ) - nhood[l].rpos );    }
               int ulen = ustop - ustart, maxid = 0, buf = 3;
               for ( int l = 0; l < ids.isize( ); l++ )
                    maxid = Max( maxid, ids[l].isize( ) );
               matrix<char> M( nhood.size( ), maxid + buf + ulen, ' ' );
               for ( int l = 0; l < nhood.isize( ); l++ )
               {    for ( int x = 0; x < ids[l].isize( ); x++ )
                         M( l, maxid - ids[l].isize( ) + x ) = ids[l][x];
                    const basevector& R 
                         = (is_frag[l] ? reads : jreads)[ nhood[l].rid ];
                    for ( int x = Max( 0, nhood[l].rpos ); x < R.isize( ); x++ )
                    {    char b = nhood[l].fw ? R[x] : 3 - R[ R.isize( ) - x - 1 ];
                         int y = nhood[l].upos + x - nhood[l].rpos - ustart;
                         if ( y >= 0 && y < ulen ) 
                              M( l, maxid + buf + y ) = as_base(b);    }    }
               cout << "\nanomaly at unibase position " << u << "." << j 
                    << " [rc unibase = " << ru << "], ";
               cout << "starting " << j - ustart << " bases to the left\n";
               for ( int x = 0; x < j - ustart + maxid + buf; x++ )
                    cout << " ";
               cout << "*" << "\n";
               M.Print(cout);    }    }

     // Compute final coverage and print.  Note that forward and reverse
     // coverage must be combined.

     cout << Date( ) << ": begin main coverage analysis" << endl;
     int coverage_print_lines = 0;
     vec< vec<int> > suspicious(nuni);
     for ( size_t u = 0; u < nuni; u++ )
     {    Bool first = True, ellipsis = False;
          int nu = unibases[u].size( ), ru = to_rc[u];
          if ( delete_close_to_better[u] )
          {    cout << "\nu = " << u << " [len=" << nu << "]"
                    << ": CLOSE TO UNIPATH WITH MUCH BETTER COVERAGE, "
                    << "RECOMMEND DELETION" << endl;
               if ( !SHOW_ALL ) continue;    }
     
          // Check for holes.

          Bool holy = False;
          for ( int j = 0; j < nu; j++ )
          {    for ( int k = 0; k < 4; k++ )
                    X[k] = calls[k][u][j] + calls[3-k][ru][nu-j-1];
               if ( Sum(X) == 0 ) holy = True;    }
          if ( holy && nu < min_to_keep_with_holes )
          {    cout << "\nu = " << u << " [len=" << nu << "]"
                    << ": HAS COVERAGE HOLE, RECOMMMEND DELETION" << endl;
               to_delete[u] = True;
               if ( !SHOW_ALL ) continue;    }

          // Check for minimum quality coverage.

          const int min_qual_cov = 10;
          const int min_size_to_waive_qual_cov = 200;
          Bool near_holy = False;
          for ( int j = 0; j < nu; j++ )
          {    for ( int k = 0; k < 4; k++ )
                    Q[k] = qcalls[k][u][j] + qcalls[3-k][ru][nu-j-1];
               if ( Max(Q) < min_qual_cov ) near_holy = True;    }
          if ( near_holy && nu < min_size_to_waive_qual_cov )
          {    cout << "\nu = " << u << " [len=" << nu << "]"
                    << ": HAS NEAR COVERAGE HOLE, RECOMMMEND DELETION" << endl;
               to_delete[u] = True;
               if ( !SHOW_ALL ) continue;    }

          // Compute coverage and print.

          for ( int j = 0; j < nu; j++ )
          {    for ( int k = 0; k < 4; k++ )
               {    X[k] = calls[k][u][j] + calls[3-k][ru][nu-j-1];
                    Q[k] = qcalls[k][u][j] + qcalls[3-k][ru][nu-j-1];    }
               int D = unibases[u][j], next_best = 0;
               for ( int l = 0; l < 4; l++ )
                    if ( l != D ) next_best = Max( next_best, Q[l] );
               if ( !SHOW_ALL && X[D] >= 2 && double(Q[D])/double(next_best) >= 3.0 )
               {    if ( !ellipsis && !first )
                    {    cout << "...\n";
                         ellipsis = True;    }
                    continue;    }
               ellipsis = False;
               if (first)
               {    cout << "\nu = " << u << " [len=" << nu << "]:\n";
                    first = False;    }
               cout << "u = " << u << ", pos = " << j << ", " << as_base(D)
                    << ", A[" << X[0] << ":" << Q[0] << "], C[" << X[1] << ":" 
                    << Q[1] << "], G[" << X[2] << ":" << Q[2]
                    << "], T[" << X[3] << ":" << Q[3] << "]";
               ++coverage_print_lines;
               if ( Sum(X) == 0 ) cout << " **********";
               if ( Q[D] != Max(Q) ) cout << " !!!!!!!!!!";
               suspicious[u].push_back(j);
               cout << endl;    }    }

     // Delete killed unibases.

     for ( size_t i = 0; i < nuni; i++ )
          if ( to_delete[i] ) unibases[i].resize(0);

     // For each unibase that was killed, assemble its reads.  Also assemble reads
     // near suspicious spots.  Also assemble reads at short terminal unipaths.

     quals.ReadAll( run_dir + "/frag_reads_filt.qualb" );
     if (USE_JUMPS)
     {    jquals.ReadAll( run_dir + "/jump_reads_filt.qualb" );
          for ( size_t i = 0; i < njreads; i++ )
               jquals[i].resize( trim_to[i] );    }

     PairsManager pairs( run_dir + "/frag_reads_filt.pairs" );
     PairsManager jpairs;
     if (USE_JUMPS)
     {    jpairs.Read( run_dir + "/jump_reads_filt.pairs" );    }

     for ( size_t u = 0; u < nuni; u++ )
     {    int nu = unibases[u].size( ), ru = to_rc[u];
          Bool little_hang = ( nexts[u].empty( ) && nu <= 200 );
          if ( !to_delete[u] && suspicious[u].empty( ) && !little_hang ) continue;
          if (PRINT_ASSEMBLIES) cout << "\n";
          cout << Date( ) << ": assembling ";
          if ( !to_delete[u] && !little_hang ) cout << "some ";
          cout << "reads from unibase " << u << endl;
          vecbasevector R;
          vecqualvector Q;
          vec<int> fwids_frag, rcids_frag, fwids_jump, rcids_jump;
          int M = K;
          for ( size_t l = SEGS_START[u]; l < SEGS_START[u+1]; l++ )
          {    int rid = SEGS[l].rid;

               Bool hits_suspicious = False;
               int upos = SEGS[l].upos;
               for ( int j = 0; j < suspicious[u].isize( ); j++ )
               {    int s = suspicious[u][j];
                    if ( s >= upos - M && s < upos + reads[rid].isize( ) + M )
                         hits_suspicious = True;    }
               if ( !to_delete[u] && !little_hang && !hits_suspicious ) continue;

               fwids_frag.push_back(rid);
	       if ( pairs.isPaired( rid ) ) 
		 rcids_frag.push_back( pairs.getPartnerID( rid ) ); }
          for ( size_t l = SEGS_START[ru]; l < SEGS_START[ru+1]; l++ )
          {    int rid = SEGS[l].rid;

               Bool hits_suspicious = False;
               int upos = SEGS[l].upos - 1;
               for ( int j = 0; j < suspicious[ru].isize( ); j++ )
               {    int s = suspicious[ru][j];
                    if ( s >= upos - M && s < upos + reads[rid].isize( ) + M )
                         hits_suspicious = True;    }
               if ( !to_delete[u] && !little_hang && !hits_suspicious ) continue;

               rcids_frag.push_back(rid);
	       if ( pairs.isPaired( rid ) ) 
		 fwids_frag.push_back( pairs.getPartnerID( rid ) ); }
          if (USE_JUMPS)
          {    for ( size_t l = JSEGS_START[u]; l < JSEGS_START[u+1]; l++ )
               {    int rid = JSEGS[l].rid;

                    Bool hits_suspicious = False;
                    int upos = JSEGS[l].upos;
                    for ( int j = 0; j < suspicious[u].isize( ); j++ )
                    {    int s = suspicious[u][j];
                         if ( s >= upos - M && s < upos + jreads[rid].isize( ) + M )
                              hits_suspicious = True;    }
                    if ( !to_delete[u] && !little_hang && !hits_suspicious ) 
                         continue;

                    fwids_jump.push_back(rid);
	            if ( jpairs.isPaired(rid) ) 
		         rcids_jump.push_back( jpairs.getPartnerID(rid) );    }
               for ( size_t l = JSEGS_START[ru]; l < JSEGS_START[ru+1]; l++ )
               {    int rid = JSEGS[l].rid;

                    Bool hits_suspicious = False;
                    int upos = JSEGS[l].upos - 1;
                    for ( int j = 0; j < suspicious[ru].isize( ); j++ )
                    {    int s = suspicious[ru][j];
                         if ( s >= upos - M && s < upos + jreads[rid].isize( ) + M )
                              hits_suspicious = True;    }
                    if ( !to_delete[u] && !little_hang && !hits_suspicious ) 
                         continue;

                    rcids_jump.push_back(rid);
	            if ( jpairs.isPaired(rid) ) 
		         fwids_jump.push_back( jpairs.getPartnerID(rid) );    }    }
          UniqueSort(fwids_frag), UniqueSort(rcids_frag);
          if (USE_JUMPS)
          {    UniqueSort(fwids_jump), UniqueSort(rcids_jump);    }
          for ( int i = 0; i < fwids_frag.isize( ); i++ )
          {    R.push_back( reads[ fwids_frag[i] ] );
               Q.push_back( quals[ fwids_frag[i] ] );    }
          if (USE_JUMPS)
          {    for ( int i = 0; i < fwids_jump.isize( ); i++ )
               {    R.push_back( jreads[ fwids_jump[i] ] );
                    Q.push_back( jquals[ fwids_jump[i] ] );    }    }
          for ( int i = 0; i < rcids_frag.isize( ); i++ )
          {    R.push_back( reads[ rcids_frag[i] ].ReverseComplement( ) );
               Q.push_back( Reverse( quals[ rcids_frag[i] ] ) );    }
          if (USE_JUMPS)
          {    for ( int i = 0; i < rcids_jump.isize( ); i++ )
               {    R.push_back( jreads[ rcids_jump[i] ].ReverseComplement( ) );
                    Q.push_back( Reverse( jquals[ rcids_jump[i] ] ) );    }    }
          HyperBasevector h;
          int basket_verbosity = 0;
          if ( BinMember( print_assemblies_verbose_list, (int) u ) ) 
               basket_verbosity = 2;
          KmerBasket( R, Q, h, 40, 10, basket_verbosity, cout );
          if (PRINT_ASSEMBLIES) cout << h;
          for ( int j = 0; j < h.EdgeObjectCount( ); j++ )
          {    if ( h.EdgeObject(j).isize( ) >= K + 1 )
                    unibases.push_back( h.EdgeObject(j) );    }    }

     // Add good new stuff.

     // unibases.Append(good_new_stuff);
     unibases.Append(new_stuff);

     // Rebuild and write unipaths.

     vec<basevector> missing_genomic_kmers;
     int64_t non_genomic_kmers = 0, N50_unibase = 0;
     if ( WRITE || VALIDATE )
     {    cout << "\n" << Date( ) << ": rebuilding unipaths" << endl;
          vecKmerPath paths, pathsrc, unipaths;
          vec<tagged_rpint> pathsdb, unipathsdb;
          ReadsToPathsCoreY( unibases, K, paths, pathsrc, pathsdb );
          Unipath( paths, pathsrc, pathsdb, unipaths, unipathsdb );
          digraph A;
          BuildUnipathAdjacencyGraph( paths, pathsrc, pathsdb, unipaths, 
               unipathsdb, A );
          HyperKmerPath h;
          BuildUnipathAdjacencyHyperKmerPath( K, A, unipaths, h );
          KmerBaseBroker kbb( K, paths, pathsrc, pathsdb, unibases );
          unibases.clear( );
          for ( int i = 0; i < h.EdgeObjectCount( ); i++ )
               unibases.push_back( kbb.Seq( h.EdgeObject(i) ) );

          // Compute summary stats from unibases.

          if (VALIDATE)
          {    cout << Date( ) << ": computing summary stats from new unibases" 
                    << endl;
               if ( K == 80 )
               {    UnibaseSummaryStats<80>( unibases, genome, 
                         missing_genomic_kmers, non_genomic_kmers, 
                         N50_unibase );    }
               else ForceAssert( 0 == 1 );    }

          // Write new files.

          if (WRITE)
          {    cout << Date( ) << ": writing new files" << endl;
               unibases.WriteAll( run_dir + "/" + OUT_HEAD + ".unibases.k80" );
               paths.WriteAll( run_dir + "/" + OUT_HEAD + ".paths.k80" );
               pathsrc.WriteAll( run_dir + "/" + OUT_HEAD + ".paths_rc.k80" );
               unipaths.WriteAll( run_dir + "/" + OUT_HEAD + ".unipaths.k80" );
               BinaryWriter::writeFile( run_dir + "/" + OUT_HEAD + ".pathsdb.k80", pathsdb );
               BinaryWriter::writeFile( run_dir + "/" + OUT_HEAD + ".unipathsdb.k80",
                    unipathsdb );    }    }

     // Print missing genomic kmers.

     if (DUMP_MISSING_GENOMIC_KMERS)
     {    cout << "\nMISSING GENOMIC KMERS:\n";
          for ( int i = 0; i < missing_genomic_kmers.isize( ); i++ )
          {    missing_genomic_kmers[i].Print( 
                    cout, "missing_" + ToString(i) );    }    }

     // Finish up.

     Remove(temp_file), Remove(temp_file2), Remove(unifile);
     cout << "\nSUMMARY STATISTICS:\n";
     cout << "coverage = " << setprecision(4) 
          << double(total)/double( unibases.sumSizes( ) ) << "\n";
     cout << "concordance = " << PERCENT_RATIO( 3, qgood, qtotal ) << "\n"
          << "deleted unipaths = " << Sum(to_delete) << "\n"
          << "coverage print lines = " << coverage_print_lines << "\n";
     if (VALIDATE)
     {    cout << "missing genomic kmers = " << missing_genomic_kmers.size( ) 
               << " [was " << old_missing_genomic_kmers.size( ) << "]" << endl;
          cout << "non-genomic kmers = " << non_genomic_kmers 
               << " [was " << old_non_genomic_kmers << "]" << endl;
          cout << "N50 unibase = " << N50_unibase 
               << " [was " << old_N50_unibase << "]" << endl;    }
     cout << "time used = " << TimeSince(clock) << "\n"
          << "\n" << Date( ) << ": done" << endl;    }
