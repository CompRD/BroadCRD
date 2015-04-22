///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// GapDumpster.  Estimate gaps by simulation.
//
// KNOWN ISSUES
//
// 1. Use of fragment pairs.  Theoretically one could increase gap accuracy by using
// fragment pairs in addition to jump pairs.  
// (a) The primary issue with this is that use of fragment pairs tends to lead to 
// incorrect gap estimates because of incorrect read placement at tandem repeats.  
// One could try to mitigate this behavior by removing placements that have better 
// placements on the unipath graph, or which land on unipaths of copy number greater 
// than ploidy.
// (b) There are also issues with the way we combine read pairs from different
// libraries, and these could be exacerbated when fragment and jump libraries are
// combined.
// (c) The variable min_sum is probably not set to the right value for fragment
// pairs.  Generally, the algorithm would need to be adapted to take account of the
// fact that fragment pairs frequently overlap.
//
// 2. Searching for the best gap value.  In a substantial fraction of cases, the
// the existing algorithm does not yield the minimum.  One can do better by more
// carefully seeking out every potential local minimum, but run time is increased
// and it is not clear that such an approach would substantially improve results.
//
// 3. Short reads.  This code ignores reads of length < 20.  This seems OK.
//
// 4. Computational performance.  This code has not been tested on large genomes
// and is unlikely to work efficiently on them without changes.
//
// 5. Evaluation.  The built-in evaluation only defines truth values for gaps that
// are flanked by unique perfect sequence.  Also in evaluation mode, the code only
// runs on gaps for which truth values have been computed.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
// MakeDepend: dependency PlotPoints

#include <omp.h>
#include <unistd.h>

#include "Basevector.h"
#include "FastIfstream.h"
#include "FeudalMimic.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "ParseSet.h"
#include "Superb.h"
#include "kmers/KmerRecord.h"
#include "math/Functions.h"
#include "paths/AssemblyCleanupTools.h"
#include "paths/PairDistCorrection.h"
#include "paths/ReadLoc.h"
#include "paths/Unipath.h"
#include "random/Random.h"
#include "system/HostName.h"

template<int K> void kmer<K>::SetToSubOf( const basevector& source,
     const size_type start )
{   size_t len = K;
    AssertLe( start, source.size() );
    AssertLe( len, source.size()-start );
    size_t end = len & ~15;
    int32_t* dst = (int32_t*) &data_;
    for ( size_t idx = 0; idx < end; idx += 16)
         *dst++ = source.extractKmer(start+idx, 16);
    if ( end < len ) *dst = source.extractKmer(start+end, len-end);   }

// template void kmer<20>::SetToSubOf( const basevector&, const size_type );

template<int K> void MakeKmerLookup( const vecbasevector& tigs,
     vec< triple<kmer<K>,int,int> >& kmers_plus )
{    vec<int64_t> starts;
     starts.push_back(0);
     for ( size_t i = 0; i < tigs.size( ); i++ )
     {    const basevector& u = tigs[i];
          starts.push_back( starts.back( ) + u.isize( ) - K + 1 );    }
     kmers_plus.resize( starts.back( ) );
     for ( size_t i = 0; i < tigs.size( ); i++ )
     {    const basevector& u = tigs[i];
          #pragma omp parallel for
          for ( int jz = 0; jz <= u.isize( ) - K; jz += 1000 )
          {    kmer<K> x, xrc;
               for ( int j = jz; j <= Min( u.isize( ) - K, jz + 1000 ); j++ )
               {    int64_t r = starts[i] + j;
                    x.SetToSubOf( u, j ); 
                    xrc = x;
                    xrc.ReverseComplement( );
                    Bool fw = ( x < xrc );
                    kmers_plus[r].first = ( fw ? x : xrc );
                    kmers_plus[r].second = i; 
                    kmers_plus[r].third = ( fw ? j : -j-1 );    }    }    }
     ParallelSort(kmers_plus);    }

double GapToKS( vec< pair<int,IntDistribution> >& computed, 
     const vec< triple<int,int,int> >& grunkle, const int gap, const vec<int>& rand,
     int& rptr, const int m1len, const int m2len, const int read_length, 
     const int max_dist, const int K, const vecbitvector& uniq, const int m1, 
     const int m2, const vec<double>& density2, const vec<int>& dists,
     const vec<int>& starts_fw, const vec<int>& stops_rc, 
     const double half_cov, const vec<int>& DISTS )
{    
     // Test to see if simulated distribution already computed.

     int did = -1;
     for ( int j = 0; j < computed.isize( ); j++ )
     {    if ( computed[j].first == gap )
          {    did = j;
               break;    }    }

     // Predict number of links.

     double exp_links = 0;
     /*
     int low1 = lower_bound( starts_fw.begin( ), starts_fw.end( ),
          m1len - max_dist ) - starts_fw.begin( );
     int high1 = starts_fw.size( );
     for ( int j = low1; j < high1; j++ )
     {    int start = starts_fw[j];
          int r = rand[rptr];
          if ( ++rptr == rand.isize( ) ) rptr = 0;
          int d = DISTS[ r % DISTS.size( ) ];
          int stop = start + d - m1len - gap;
          if ( stop < read_length || stop >= m2len ) continue;

          // Compare to expected coverage on m2.

          int delta = 100;
          int count2 = 0;
          int low = lower_bound( stops_rc.begin( ), stops_rc.end( ), stop - delta )
               - stops_rc.begin( );
          int high = lower_bound( stops_rc.begin( ), stops_rc.end( ), stop + delta )
               - stops_rc.begin( );
          for ( int k = low; k < high; k++ )
               if ( Abs( stop - stops_rc[k] ) <= delta ) count2++;
          double exp2 = half_cov * double( 2*delta + 1 );
          double add = Min( 1.0, double(count2) / exp2 );
          exp_links += add;    }
     */
     /*
     cout << "gap = " << gap << ", expected links = " // ZZZZZZZZZZZZZZZZZZZZZZZZZZZ
          << exp_links << ", found " << 2*dists.size( ) // ZZZZZZZZZZZZZZZZZZZZZZZZZ
          << endl; // ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
     */

     // Compute simulated distribution.

     vec<int> X;
     if ( did < 0 )
     {
          // Make sure that:
          // read_length <= stop2
          // and            stop2 < max_dist
          // and            stop2 <= m2len.

          int min_sum = read_length + gap + m1len;
          // int min_sum = K + gap + m1len;
          while(1)
          {    int stop2 = min_sum - m1len - gap;
               if ( stop2-K >= m2len ) break;
               if ( uniq[m2][stop2-K] ) break;
               min_sum++;    }
          int low = lower_bound( grunkle.begin( ), grunkle.end( ),
               make_triple( min_sum, 0, 0 ) ) - grunkle.begin( );
          int max_sum = Min( m2len, max_dist - 1 ) + m1len + gap;
          int high = upper_bound( grunkle.begin( ), grunkle.end( ),
               make_triple( max_sum+1, 0, 0 ) ) - grunkle.begin( );

          // Some complicated stuff follows.  Basically it's about pretesting to
          // see if we should make a list in advance of all the 'usable' positions.

          vec<int> usable;
          int ucount = 0, block = Max( 1, (high-low)/100 );
          for ( int m = low; m < high; m += block )
          {    int start1 = grunkle[m].second, x = grunkle[m].third;
               int stop2 = x - ( m1len - start1 ) - gap;
               if ( !uniq[m2][ stop2 - K ] ) continue;
               if ( density2[stop2] == 0 ) continue;
               ucount++;    }
          if ( ucount < 20 )
          {    for ( int m = low; m < high; m++ )
               {    int start1 = grunkle[m].second, x = grunkle[m].third;
                    int stop2 = x - ( m1len - start1 ) - gap;
                    if ( !uniq[m2][ stop2 - K ] ) continue;
                    if ( density2[stop2] == 0 ) continue;
                    usable.push_back(m);    }
               if ( usable.empty( ) ) return 2;    }

          // Set up for main loop.

          const int sample = 20000;
          int iteration = 0;
          while(1)
          {    
               // Pick start1 and x.

               int m;
               if ( usable.empty( ) )
               {    if ( high - low <= sample )
                    {    m = low + iteration++;
                         if ( m == high ) break;    }
                    else
                    {    int r = rand[rptr];
                         if ( ++rptr == rand.isize( ) ) rptr = 0;
                         m = low + ( r % (high-low) );    }    }
               else if ( usable.isize( ) <= sample )
               {    if ( iteration == usable.isize( ) ) break;
                    m = usable[iteration++];    }
               else
               {    int r = rand[rptr];
                    if ( ++rptr == rand.isize( ) ) rptr = 0;
                    m = usable[ r % usable.size( ) ];    }
               int start1 = grunkle[m].second, x = grunkle[m].third;
     
               // Find stop2 and test.
     
               int stop2 = x - ( m1len - start1 ) - gap;
               if ( !uniq[m2][ stop2 - K ] ) continue;
               if ( density2[stop2] == 0 ) continue;
               if ( density2[stop2] < 1 )
               {    int r = rand[rptr];
                    if ( ++rptr == rand.isize( ) ) rptr = 0;
                    if ( double( r % 100 ) / 100.0 > density2[stop2] ) continue;    }

               // Accept x, smooth the distribution.

               for ( int j = 0; j <= 20; j++ )
                    X.push_back( x + j - 10 );
               if ( X.isize( ) == sample*21 ) break;    }
          if ( X.empty( ) ) return 2;    }

     // Define two distributions and compare them.

     IntDistribution dg;
     if ( did < 0 ) 
     {    dg.from_hits(X);
          computed.push( gap, dg );    }
     else dg = computed[did].second;

     vec<int> distsb(dists);
     int dtop;
     for ( dtop = 0; dtop < dists.isize( ); dtop++ )
          if ( dists[dtop] + gap >= max_dist * 2 ) break;
     distsb.resize(dtop);
     if ( distsb.empty( ) ) return 2;

     double KS = 0.0; 
     int dpos = 0;
     // for ( int u = 0; u < max_dist; u++ )
     for ( int u = 0; u < max_dist * 2; u++ )
     {    while( dpos < distsb.isize( ) && distsb[dpos] + gap <= u ) dpos++;
          double pd = double(dpos) / double( distsb.size( ) );
          // KS = Max( KS, Abs( dg.prob_le(u) - pd ) );    }
          KS += Abs( dg.prob_le(u) - pd );    }
     KS /= double(max_dist);
     /*
     KS *= 1.0 - 0.02 * 
           log2( Abs( exp_links - 2*dists.isize( ) ) 
                / Max( exp_links, double( 2*dists.isize( ) ) ) );
     */
     return KS;    }

void PrintSideBySide( const String& z1, const String& z2, const int N )
{    istringstream zin1(z1), zin2(z2);
     String line1, line2;
     while(1)
     {    getline( zin1, line1 );
          getline( zin2, line2 );
          if ( !zin1.fail( ) && !zin2.fail( ) )
          {    cout << line1;
               for ( int i = 0; i < (int) N - (int) line1.size( ); i++ )
                    cout << " ";
               cout << line2 << "\n";    }
          else if ( zin1.fail( ) && !zin2.fail( ) )
          {    for ( int i = 0; i < N; i++ )
                    cout << " ";
               cout << line2 << "\n";    }
          else break;    }    }

template<int K> void GetBounds( const vec< kmer<K> >& kmers, kmer<K>& x,
     kmer<K>& xrc, Bool& fw, int64_t& low, int64_t& high )
{    xrc = x;
     xrc.ReverseComplement( );
     fw = ( x < xrc );
     if ( !fw ) x = xrc;
     low = lower_bound( kmers.begin( ), kmers.end( ), x ) - kmers.begin( );
     high = upper_bound( kmers.begin( ), kmers.end( ), x ) - kmers.begin( );    }

// GetStartStop.  Find start-stop points of reads.  The algorithm assumes that jumps 
// have been reverse complemented, which is not something we should rely on.

template<int K, int K2> void GetStartStop( const String& run_dir, 
     const String& head, const int ntigs, 
     const vec< triple<kmer<K>,int,int> >& kmers_plus, const vec< kmer<K> >& kmers, 
     const vec< triple<kmer<K2>,int,int> >& kmers_plus2,
     const vec< kmer<K2> >& kmers2, vec< vec< vec<int> > >& starts_fw, 
     vec< vec< vec<int> > >& stops_rc, const int VERBOSITY, vec<int>& read_length )
{
     cout << Date( ) << ": loading " << head << endl;
     String pairs_file;
     if ( head == "frag_reads_filt" ) 
          pairs_file = run_dir + "/" + head + "_cpd.pairs";
     else pairs_file = run_dir + "/" + head + ".pairs";
     PairsManager pairs(pairs_file);
     pairs.makeCache( );
     vecbasevector bases( run_dir + "/" + head + ".fastb" );
     int nlibs = starts_fw.size( ), nlibs_this = pairs.nLibraries( );
     if ( VERBOSITY >= 1 )
     {    for ( int l = 0; l < nlibs_this; l++ )
          {    cout << Date( ) << ": library " << nlibs + l << " is of type "
                    << head.Before( "_filt" ) << "\n";    }    }
     starts_fw.resize( nlibs + nlibs_this ), stops_rc.resize( nlibs + nlibs_this );
     for ( int l = 0; l < nlibs_this; l++ )
          read_length.push_back(0);
     for ( size_t id = 0; id < bases.size( ); id++ )
     {    int64_t pid = pairs.getPairID(id);
          int lib = pairs.libraryID(pid);
          read_length[nlibs+lib] = Max( read_length[nlibs+lib], bases[id].isize( ) );
               }
     for ( int l = 0; l < nlibs_this; l++ )
     {    starts_fw[nlibs+l].resize(ntigs), stops_rc[nlibs+l].resize(ntigs);    }
     cout << Date( ) << ": looking up " << head << endl;
     const size_t batch_size = 100000;
     #pragma omp parallel for
     for ( size_t xx = 0; xx < bases.size( ); xx += batch_size )
     {    kmer<K> x, xrc;
          kmer<K2> x2, x2rc;
          int count = 0;
          int64_t low, high;
          Bool fw, fw_tig;
          int tig, pos;
          vec< triple<int,int,int> > starts_fw_pushes, stops_rc_pushes;
          for ( size_t id = xx; id < Min( bases.size( ), xx + batch_size ); id++ )
          {    int64_t pid = pairs.getPairID(id);
               int lib = pairs.libraryID(pid);

               // Ignore tiny reads.

               if ( bases[id].isize( ) < K2 ) continue;

               // First process the short reads and then the longer reads.  Note
               // that for the longer jumps, we take the last kmer of the reads,
               // under the assumption that they have been reverse complemented.
               // (This is horrible.)  Also, if there is more than one placement, 
               // we make a pseudo-random choice.

               if ( bases[id].isize( ) < K )
               {    x2.SetToSubOf( bases[id], 0 );
                    GetBounds( kmers2, x2, x2rc, fw, low, high );
                    if ( high - low >= 1 )
                    {    int j = low + ( (count++) % (high-low) );
                         tig = kmers_plus2[j].second, pos = kmers_plus2[j].third;
                         fw_tig = ( pos >= 0 );
                         if ( !fw_tig ) pos = -pos-1;
                         if ( !fw_tig ) pos += K2;
                         if ( fw == fw_tig ) starts_fw_pushes.push( lib, tig, pos );
                         else stops_rc_pushes.push( lib, tig, pos );    }    }
               else
               {    if ( head.Contains( "jump", 0 ) )
                         x.SetToSubOf( bases[id], bases[id].isize( ) - K );
                    else x.SetToSubOf( bases[id], 0 );
                    GetBounds( kmers, x, xrc, fw, low, high );
                    if ( high - low >= 1 )
                    {    int j = low + ( (count++) % (high-low) );
                         tig = kmers_plus[j].second, pos = kmers_plus[j].third;
                         fw_tig = ( pos >= 0 );
                         if ( !fw_tig ) pos = -pos-1;
                         if ( head.Contains( "jump", 0 ) )
                         {    if (fw_tig) pos -= ( bases[id].isize( ) - K );
                              else pos += bases[id].isize( );     }
                         else
                         {    if ( !fw_tig ) pos += K;    }
                         if ( fw == fw_tig ) starts_fw_pushes.push( lib, tig, pos );
                         else stops_rc_pushes.push( lib, tig, pos );    }    }    }
          #pragma omp critical
          {    for ( int j = 0; j < starts_fw_pushes.isize( ); j++ )
               {    const triple<int,int,int>& x = starts_fw_pushes[j];
                    starts_fw[nlibs+x.first][x.second].push_back(x.third);    }
               for ( int j = 0; j < stops_rc_pushes.isize( ); j++ )
               {    const triple<int,int,int>& x = stops_rc_pushes[j];
                    stops_rc[nlibs+x.first][x.second].push_back(x.third);   }   }   }
     cout << Date( ) << ": sorting" << endl;
     for ( int l = 0; l < nlibs_this; l++ )
     {    for ( int t = 0; t < ntigs; t++ )
          {    Sort( starts_fw[nlibs+l][t] ); 
               Sort( stops_rc[nlibs+l][t] );    }    }    }

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(SUBDIR, "test");
     CommandArgument_String(SCAFFOLDS_IN);
     CommandArgument_Int_OrDefault_Doc(NUM_THREADS, -1,
          "number of threads to be used (use all available if negative)");
     CommandArgument_String_OrDefault_Doc(TIGS, "",
          "if specified, only process gaps have these contigs on their left");
     CommandArgument_Int_OrDefault_Doc(VERBOSITY, -1,
          "VERBOSITY=1: print one line for each gap\n"
          "VERBOSITY=2: print one line for each gap estimate\n"
          "VERBOSITY=3: print pair placements\n"
          "(Note VERBOSITY>=2 should only be used if TIGS is set to one contig.)");
     CommandArgument_Bool_OrDefault(WRITE, True);
     CommandArgument_Bool_OrDefault(VALIDATE, True);
     CommandArgument_Bool_OrDefault(SHOW_START_STOP, False);
     CommandArgument_String_OrDefault_Doc(PLOT_HEAD, "",
          "if specified, generate PLOT_HEAD.png" );
     CommandArgument_String_OrDefault_Doc(LIB_TYPES, "jump",
          "theoretically could be {frag,jump,long}, see comments in code");
     CommandArgument_String_OrDefault(LIBS, "");
     EndCommandArguments;

     // Define directories, etc.

     cout << Date( ) << ": begin" << endl;
     double clock = WallClockTime( );
     String data_dir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
     String head = sub_dir + "/" + SCAFFOLDS_IN;
     vec<String> lib_types_to_use;
     ParseStringSet( LIB_TYPES, lib_types_to_use );
     vec<int> libs_to_use;
     ParseIntSet( LIBS, libs_to_use );

     // Thread control

     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );

     // Load assembly.

     String supers_file = sub_dir + "/" + SCAFFOLDS_IN + ".superb";
     vecbasevector tigs( sub_dir + "/" + SCAFFOLDS_IN + ".contigs.fastb" );
     int ntigs = tigs.size( );
     vec<superb> scaffolds;
     ReadSuperbs( supers_file, scaffolds );
     const int ploidy = FirstLineOfFile( data_dir + "/ploidy" ).Int( );
     read_locs_on_disk locs_file( head, run_dir );
     vec<int> tigs_to_process;
     if ( TIGS != "" ) 
     {    ParseIntSet( TIGS, tigs_to_process );
          if ( VERBOSITY < 0 ) VERBOSITY = 1;    }
     if ( VERBOSITY >= 2 && !tigs_to_process.solo( ) )
     {    cout << "If you use VERBOSITY >= 2, then you need to specify "
               << "a single contig.\n";
          exit(1);    }

     // Load genome.

     vecbasevector genome;
     if (VALIDATE) genome.ReadAll( data_dir + "/genome.fastb" );

     // Find true gaps.

     vec<int> true_gap;
     vec<Bool> true_gap_computed;
     if (VALIDATE)
     {    cout << Date( ) << ": computing gkmers_plus" << endl;
          const int L = 1000;
          vec< triple<kmer<L>,int,int> > gkmers_plus;
          MakeKmerLookup( genome, gkmers_plus );
          cout << Date( ) << ": forming gkmers" << endl;
          vec< kmer<L> > gkmers( gkmers_plus.size( ) );
          for ( size_t i = 0; i < gkmers.size( ); i++ )
               gkmers[i] = gkmers_plus[i].first;
          cout << Date( ) << ": finding true gaps" << endl;
          for ( int s = 0; s < scaffolds.isize( ); s++ )
          {    const superb& S = scaffolds[s];
               for ( int j = 0; j < S.Ngaps( ); j++ )
               {    int m1 = S.Tig(j), m2 = S.Tig(j+1);
                    if ( ( TIGS != "" && !BinMember( tigs_to_process, m1 ) )
                         || tigs[m1].isize( ) < L || tigs[m2].isize( ) < L )
                    {    true_gap.push_back(0);
                         true_gap_computed.push_back(False);
                         continue;    }
                    kmer<L> x, xrc;
                    x.SetToSubOf( tigs[m1], tigs[m1].isize( ) - L );
                    xrc = x;
                    xrc.ReverseComplement( );
                    Bool fw1 = ( x < xrc );
                    if ( !fw1 ) x = xrc;
                    int64_t low1 = lower_bound( gkmers.begin( ), gkmers.end( ), x )
                         - gkmers.begin( );
                    int64_t high1 = upper_bound( gkmers.begin( ), gkmers.end( ), x )
                         - gkmers.begin( );
                    if ( high1 != low1 + 1 )
                    {    true_gap.push_back(0);
                         true_gap_computed.push_back(False);
                         continue;    }
                    x.SetToSubOf( tigs[m2], 0 );
                    xrc = x;
                    xrc.ReverseComplement( );
                    Bool fw2 = ( x < xrc );
                    if ( !fw2 ) x = xrc;
                    int64_t low2 = lower_bound( gkmers.begin( ), gkmers.end( ), x )
                         - gkmers.begin( );
                    int64_t high2 = upper_bound( gkmers.begin( ), gkmers.end( ), x )
                         - gkmers.begin( );
                    if ( high2 - low2 != 1
                         || gkmers_plus[low1].second != gkmers_plus[low2].second )
                    {    true_gap.push_back(0);
                         true_gap_computed.push_back(False);
                         continue;    }
                    int pos1 = gkmers_plus[low1].third; 
                    int pos2 = gkmers_plus[low2].third;
                    Bool xfw1 = ( pos1 >= 0 ), xfw2 = ( pos2 >= 0 );
                    if ( !xfw1 ) pos1 = -pos1-1;
                    if ( !xfw2 ) pos2 = -pos2-1;
                    const int true_gap_upper_bound = 100000;
                    int g = ( fw1 == xfw1 ? pos2 - (pos1+L) : pos1 - (pos2+L) );
                    if ( ( int(fw1) + int(fw2) + int(xfw1) + int(xfw2) ) % 2 != 0
                         || g > true_gap_upper_bound )
                    {    // PRINT8( m1, m2, int(fw1), int(fw2), 
                         // int(xfw1), int(xfw2), pos1, pos2 );
                         true_gap.push_back(0);
                         true_gap_computed.push_back(False);
                         continue;    }
                    // cout << "gap from " << m1 << " to " 
                    //      << m2 << " is " << g << endl;
                    true_gap.push_back(g);
                    true_gap_computed.push_back(True);    }    }
          cout << Date( ) << ": found true size for " << Sum(true_gap_computed)
               << " gaps of " << true_gap_computed.size( ) << endl;    }

     // Build a map of kmers in the contigs.

     cout << Date( ) << ": building kmer maps for contigs" << endl;
     const int K = 40;
     const int K2 = 20;
     vec< triple<kmer<K>,int,int> > kmers_plus;
     MakeKmerLookup( tigs, kmers_plus );
     vec< triple<kmer<K2>,int,int> > kmers_plus2;
     MakeKmerLookup( tigs, kmers_plus2 );

     // Find unique kmers in the contigs.  

     cout << Date( ) << ": finding unique kmers" << endl;
     vecbitvector uniq;
     {    Mimic( tigs, uniq );
          for ( size_t i = 0; i < kmers_plus.size( ); i++ )
          {    size_t j;
               for ( j = i + 1; j < kmers_plus.size( ); j++ )
                    if ( kmers_plus[j].first != kmers_plus[i].first ) break;
               if ( j - i == 1 )
               {    int pos = kmers_plus[i].third;
                    if ( pos < 0 ) pos = -pos-1;
                    uniq[ kmers_plus[i].second ].Set( pos, True );    }
               i = j - 1;    }    }

     // Form kmers.

     cout << Date( ) << ": copying" << endl;
     vec< kmer<K> > kmers( kmers_plus.size( ) );
     for ( size_t i = 0; i < kmers.size( ); i++ )
          kmers[i] = kmers_plus[i].first;
     vec< kmer<K2> > kmers2( kmers_plus2.size( ) );
     for ( size_t i = 0; i < kmers2.size( ); i++ )
          kmers2[i] = kmers_plus2[i].first;

     // Find oriented start points in "the" jumping library.  
     // WARNING: THIS ASSUMES THAT READS HAVE BEEN
     // REVERSE COMPLEMENTED.  WARNING: THIS ASSUMES THAT READS HAVE SIZE AT LEAST K.
     // WARNING: I'M NOT SURE THAT THE FW/RC BOOKKEEPING HERE IS RIGHT.
     // NOTE WE'RE NOT USING LONG JUMPS.
     //
     // Also get read lengths, actually the maximum read length for each library.
     // It isn't really right to use one value for each library.

     vec< vec< vec<int> > > starts_fw, stops_rc;
     vec<int> read_length;
     if ( Member( lib_types_to_use, String("frag") ) )
     {    GetStartStop( run_dir, "frag_reads_filt", ntigs, kmers_plus, kmers, 
               kmers_plus2, kmers2, starts_fw, stops_rc, VERBOSITY,
               read_length );    }
     int nlibs_frag = starts_fw.size( );
     if ( Member( lib_types_to_use, String("jump") ) )
     {    GetStartStop( run_dir, "jump_reads_filt", ntigs, kmers_plus, kmers, 
               kmers_plus2, kmers2, starts_fw, stops_rc, VERBOSITY,
               read_length );    }
     int nlibs_jump = starts_fw.size( ) - nlibs_frag;
     if ( Member( lib_types_to_use, String("long") ) )
     {    if ( IsRegularFile( run_dir + "/long_jump_reads_filt.fastb" ) )
          {    GetStartStop( run_dir, "long_jump_reads_filt", ntigs, kmers_plus,
                    kmers, kmers_plus2, kmers2, starts_fw, stops_rc, VERBOSITY,
                    read_length );    }     }
     int nlibs_long = starts_fw.size( ) - nlibs_frag - nlibs_jump;
     int nlibs = starts_fw.size( );

     // Compute distributions.

     vec<IntDistribution> jdist;
     /*
     {    String distribs_file = run_dir + "/jump_reads_filt_cpd.distribs";
          vec<IntDistribution> distribs;
          BinaryReader::readFile( distribs_file.c_str( ), &distribs );
          for ( size_t i = 0; i < distribs.size( ); i++ ) 
          {    IntDistribution dist_lt, dist_gt;
               distribs[i].split( &dist_lt,  &dist_gt );
               jdist.push_back(dist_gt);    }    }
     cout << Date( ) << ": ";
     PRINT( jdist.size( ) );
     */

     // Set a lower bound for contigs to use.  It's not clear that this makes sense.

     vec<int> len;
     int64_t total = 0, part = 0, min_to_use = 0;
     for ( size_t m = 0; m < tigs.size( ); m++ )
     {    len.push_back( tigs[m].size( ) );
          total += tigs[m].size( );    }
     ReverseSort(len);
     for ( size_t i = 0; i < len.size( ); i++ )
     {    part += len[i];
          min_to_use = len[i];
          if ( double(part)/double(total) >= 0.5 ) break;    }
     cout << Date( ) << ": "; PRINT(min_to_use);

     // Define visibility at given distances.

     int max_tig = 0;
     for ( size_t m = 0; m < tigs.size( ); m++ )
          max_tig = Max( max_tig, tigs[m].isize( ) );
     cout << Date( ) << ": maximum contig size = "
          << ToStringAddCommas(max_tig) << endl;
     vec<int64_t> VISIBLE( max_tig, 0 );
     for ( size_t m = 0; m < tigs.size( ); m++ )
     {    for ( size_t j = 0; j < tigs[m].size( ); j++ )
               VISIBLE[j] += tigs[m].isize( ) - j;    }

     // Compute distributions for libraries and upper cutoffs for them.  There are
     // two passes: first we find the cutoffs, then we compute the distributions,
     // taking into account the cutoffs.

     cout << Date( ) << ": computing library distributions" << endl;
     vec<Bool> lfail( nlibs, False );
     vec<int> max_dist(nlibs);
     vec< vec<int> > DISTS(nlibs);
     for ( int pass = 1; pass <= 2; pass++ )
     {    for ( int l = 0; l < nlibs; l++ )
               DISTS[l].clear( );
          for ( size_t m = 0; m < tigs.size( ); m++ )
          {    if ( tigs[m].isize( ) < min_to_use ) continue;
               int n = tigs[m].size( );
               vec<read_loc> locs;
               locs_file.LoadContig( m, locs );
               for ( int i = 0; i < locs.isize( ); i++ )
               {    const read_loc& rl = locs[i];
                    if ( !rl.PartnerPlaced( ) ) continue;
                    if ( rl.PartnerContigId( ) != m ) continue;
                    if ( !rl.Fw( ) || !rl.PartnerRc( ) ) continue;
                    if ( rl.Start( ) < 0 || rl.PartnerStop( ) > n ) continue;
                    int dist = rl.PartnerStop( ) - rl.Start( );
                    if ( dist < 0 || dist >= max_tig ) continue;
                    if ( rl.Frag( ) && !Member( lib_types_to_use, String("frag") ) ) 
                         continue;
                    if ( rl.Jump( ) && !Member( lib_types_to_use, String("jump") ) ) 
                         continue;
                    if ( rl.LongJump( ) 
                         && !Member( lib_types_to_use, String("long") ) ) 
                    {    continue;    }
                    int libid = rl.LibraryId( );
                    if ( rl.Jump( ) || rl.LongJump( ) ) libid += nlibs_frag;
                    if ( rl.LongJump( ) ) libid += nlibs_jump;
                    if ( pass == 2 && n - rl.Start( ) < max_dist[libid] ) continue;
                    if ( pass == 2 && dist >= max_dist[libid] ) continue;
                    DISTS[libid].push_back(dist);    }    }
          for ( int l = 0; l < nlibs; l++ )
               Sort( DISTS[l] );
          if ( pass == 2 ) break;

          // Define upper cutoff.  To do this, walk backwards from the end of the 
          // distribution, sampling 100 points at a time.  Taking into account
          // the visibility, estimate the noise rate and error bars for it, assuming
          // that the 100 points are all noise.  Now keep walking until we think
          // the signal:noise ratio is at least 10.  This defines the upper cutff.

          const int bin_count = 100;
          const double err_mult = 2.0;
          const double max_noise_frac = 0.1;
          for ( int l = 0; l < nlibs; l++ )
          {    if ( DISTS[l].isize( ) < bin_count ) 
               {    lfail[l] = True;
                    continue;    }
               int N = DISTS[l].size( );
               double min_freq_high = 1000000000;
               lfail[l] = True;
               for ( int rx = 0; N - (rx+1)*bin_count >= 0; rx++ )
               {    double freq = 0, effective_size = 0;
                    int high = DISTS[l][ N - (rx)*bin_count - 1 ];
                    for ( int j = N - (rx+1)*bin_count; j < N - (rx)*bin_count; j++ )
                    {    int p = DISTS[l][j], q = high;
                         freq += Mean(VISIBLE) / double( VISIBLE[p] );
                         effective_size += 
                              double(VISIBLE[q]) / double( VISIBLE[p] );    }
                    int low = DISTS[l][ N - (rx+1)*bin_count ];
                    freq /= double( high - low );
                    double freq_err = freq / sqrt(effective_size);
                    double freq_low = Max( 0.0, freq - err_mult * freq_err );
                    double freq_high = freq + err_mult * freq_err;
                    min_freq_high = Min( min_freq_high, freq_high );
                    double delta = freq_low / min_freq_high;
                    // PRINT5( low, high, freq_low, freq_high, delta );
                    if ( 1.0 / delta <= max_noise_frac )
                    {    max_dist[l] = high;
                         lfail[l] = False;
                         break;    }    }    }    }
     vec<double> lib_mean(nlibs), lib_dev(nlibs);
     for ( int l = 0; l < nlibs; l++ )
     {    if ( !lfail[l] )
          {    lib_mean[l] = Mean( DISTS[l] );
               lib_dev[l] = StdDev( DISTS[l], lib_mean[l] );    }
          cout << Date( ) << ": library " << l << ", ";
          if ( lfail[l] ) cout << "can't find upper bound" << endl;
          else
          {    cout << lib_mean[l] << " +/- " << lib_dev[l] 
                    << ", upper bound = " << max_dist[l] << "\n";    }    }

     // Compute coverage.  This does not handle innies properly.

     vec<int64_t> fw_placed( nlibs, 0 ), assayed( nlibs, 0 );
     cout << Date( ) << ": computing coverage" << endl;
     for ( size_t m = 0; m < tigs.size( ); m++ )
     {    vec<read_loc> locs;
          locs_file.LoadContig( m, locs );
          for ( int l = 0; l < nlibs; l++ )
          {    int n = tigs[m].size( );
               if ( n < max_dist[l] ) continue;
               assayed[l] += n - max_dist[l];
               for ( int i = 0; i < locs.isize( ); i++ )
               {    const read_loc& rl = locs[i];
                    if ( !rl.Fw( ) ) continue;
                    if ( !rl.PartnerPlaced( ) ) continue;
                    if ( rl.PartnerContigId( ) != m ) continue;
                    if ( !rl.Fw( ) || !rl.PartnerRc( ) ) continue;
                    if ( rl.Start( ) < 0 || rl.PartnerStop( ) > n ) continue;
                    if ( rl.Frag( ) && !Member( lib_types_to_use, String("frag") ) ) 
                         continue;
                    if ( rl.Jump( ) && !Member( lib_types_to_use, String("jump") ) ) 
                         continue;
                    if ( rl.LongJump( ) 
                         && !Member( lib_types_to_use, String("long") ) ) 
                    {    continue;    }
                    int libid = rl.LibraryId( );
                    if ( rl.Jump( ) || rl.LongJump( ) ) libid += nlibs_frag;
                    if ( rl.LongJump( ) ) libid += nlibs_jump;
                    if ( libid != l ) continue;
                    if ( n - rl.Start( ) < max_dist[l] ) continue;
                    int dist = rl.PartnerStop( ) - rl.Start( );
                    if ( dist < 0 || dist > max_dist[l] ) continue;
                    fw_placed[l]++;    }    }    }
     vec<double> half_cov(nlibs);
     for ( int l = 0; l < nlibs; l++ )
          half_cov[l] = double(fw_placed[l]) / double(assayed[l]);

     // Precompute random numbers.

     cout << Date( ) << ": precomputing random numbers" << endl;
     vec<int> rand;
     for ( int j = 0; j < 1000000; j++ )
          rand.push_back( randomx( ) );

     // Precompute values from distribution.

     cout << Date( ) << ": precomputing values" << endl;
     vec< vec<int> > D2vals(nlibs);
     for ( int l = 0; l < nlibs; l++ )
     {    if ( lfail[l] ) continue;
          int rptr = 0;
          for ( int j = 0; j < 200000; j++ )
          {    rptr++;
               if ( rptr == rand.isize( ) ) rptr = 0;
               D2vals[l].push_back( 
                    DISTS[l][ rand[rptr] % DISTS[l].size( ) ] );    }    }

     // Go through the scaffolds.

     cout << Date( ) << ": going through scaffolds" << endl;
     if ( VERBOSITY >= 1 ) cout << "\n";
     vec< triple<int,int,int> > results;
     vec<int> devs, adevs;
     vec<double> offbys, aoffbys;
     vec< pair<int,int> > gaps;
     for ( int s = 0; s < scaffolds.isize( ); s++ )
     {    const superb& S = scaffolds[s];
          for ( int jg = 0; jg < S.Ngaps( ); jg++ )
               gaps.push( s, jg );    }
     int fails = 0;
     #pragma omp parallel for schedule( dynamic, 1 )
     for ( int gi = 0; gi < gaps.isize( ); gi++ )
     {    int s = gaps[gi].first, jg = gaps[gi].second;
          const superb& S = scaffolds[s];
          int m1 = S.Tig(jg), m2 = S.Tig(jg+1);
          if ( TIGS != "" && !BinMember( tigs_to_process, m1 ) ) continue;
          int tgap = 0;
          if (VALIDATE)
          {    if ( !true_gap_computed[gi] ) continue;
               tgap = true_gap[gi];    }
          if (SHOW_START_STOP)
          {    
               #pragma omp critical
               {    cout << "start " << m1 << endl;    }    }

          // Test for plausibility of negative gaps.  Accept overlap o if at least
          // (o-K+1)/2 kmer hits are found within 10% of o.  A danger of this is 
          // that if a contig has a hanging end which "overlaps" the other contig,
          // we don't report a negative overlap.  Another problem is that this
          // won't see contigs that are out of order.
          //
          // More serious problem: gaps near the floor have compressed standard
          // deviation.

          const double req_over_frac = 0.5;
          const double over_delta = 0.1;
          const int overlap_fudge = 40;
          vec<int> overlaps, accepted_overlaps;
          for ( int p2 = 0; p2 <= Min( Max(max_dist), tigs[m2].isize( ) ) - K; p2++ )
          {    kmer<K> x, xrc;
               x.SetToSubOf( tigs[m2], p2 );
               xrc = x;
               xrc.ReverseComplement( );
               Bool fw2 = ( x < xrc );
               if ( !fw2 ) x = xrc;
               int64_t low = lower_bound( kmers.begin( ), kmers.end( ), x )
                    - kmers.begin( );
               int64_t high = upper_bound( kmers.begin( ), kmers.end( ), x )
                    - kmers.begin( );
               for ( int64_t l = low; l < high; l++ )
               {    int m = kmers_plus[l].second;
                    if ( m != m1 ) continue;
                    int p1 = kmers_plus[l].third;
                    Bool fw1 = ( p1 >= 0 );
                    if ( fw1 != fw2 ) continue;
                    if ( !fw1 ) p1 = -p1-1;
                    overlaps.push_back( p2 + ( tigs[m1].isize( ) - p1 ) );    }    }
          Sort(overlaps);
          for ( int j = 0; j < overlaps.isize( ); j++ )
          {    int o = overlaps[j];
               int support = 1;
               for ( int l = j - 1; l >= 0; l-- )
               {    if ( o - overlaps[l] > double(o) * over_delta ) break;
                    support++;    }
               for ( int l = j + 1; l < overlaps.isize( ); l++ )
               {    if ( overlaps[l] - o > double(o) * over_delta ) break;
                    support++;    }
               if ( support >= (o-K+1) * req_over_frac )
                    accepted_overlaps.push_back(o);
               int k = overlaps.NextDiff(j);
               j = k - 1;    }
          int max_overlap = overlap_fudge;
          if ( accepted_overlaps.nonempty( ) )
               max_overlap += accepted_overlaps.back( );
          if ( VERBOSITY >= 2 ) PRINT(max_overlap);

          // Find links from m1 to m2.

          vec<read_loc> locs;
          #pragma omp critical
          {    locs_file.LoadContig( m1, locs );    }
          vec< vec<int> > dists(nlibs);
          ofstream* plout = 0;
          if ( PLOT_HEAD != "" ) 
               plout = new ofstream( ( PLOT_HEAD + ".txt" ).c_str( ) );
          for ( int i = 0; i < locs.isize( ); i++ )
          {    const read_loc& rl = locs[i];
               if ( !rl.PartnerPlaced( ) ) continue;
               if ( (int) rl.PartnerContigId( ) != m2 ) continue;
               if ( !rl.Fw( ) || !rl.PartnerRc( ) ) continue;
               if ( rl.Start( ) < 0 || rl.PartnerStop( ) > tigs[m2].isize( ) )
                    continue;
               int dist1 = tigs[m1].isize( ) - rl.Start( );
               int dist2 = rl.PartnerStop( );
               int dist = dist1 + dist2;
               if ( rl.Frag( ) && !Member( lib_types_to_use, String("frag") ) ) 
                    continue;
               if ( rl.Jump( ) && !Member( lib_types_to_use, String("jump") ) ) 
                    continue;
               if ( rl.LongJump( ) && !Member( lib_types_to_use, String("long") ) ) 
                    continue;
               int libid = rl.LibraryId( );
               if ( libs_to_use.nonempty( ) && !Member( libs_to_use, libid ) )
                    continue;
               if ( rl.Jump( ) || rl.LongJump( ) ) libid += nlibs_frag;
               if ( rl.LongJump( ) ) libid += nlibs_jump;
               if ( lfail[libid] || dist < 0 ) continue;
               /*
               if ( dist >= int( ceil( double(max_dist[libid]) * 1.2 ) )  
                    + max_overlap )
               {    continue;    }
               */
               if ( PLOT_HEAD != "" ) *plout << dist1 << " " << dist2 << "\n";
               if ( VERBOSITY >= 3 ) PRINT4( libid, dist1, dist2, dist1+dist2 );
               dists[libid].push_back(dist);    }
          if ( PLOT_HEAD != "" )
          {    plout->close( );
               SystemSucceed( "PlotPoints IN=" + PLOT_HEAD + ".txt OUT="
                    + PLOT_HEAD + ".png X_AXIS_EXTEND=0 Y_AXIS_EXTEND=0 "
                    + "X_AXIS_OFFSET=0 PR=1 XMIN=0 MIN_Y=0 " 
                    + "TITLE=\"Gap distances for " + ToString(m1) + " --> "
                         + ToString(m2) + "\"" );    }

          // Give up if all libraries are essentially empty.  Remove low-power
          // libraries.

          Bool have_data = False;
          vec<double> power(nlibs);
          for ( int l = 0; l < nlibs; l++ )
          {    if ( lfail[l] ) continue;
               Sort( dists[l] );
               power[l] = lib_dev[l] / sqrt( double( dists[l].size( ) ) );
               if ( VERBOSITY >= 2 ) 
               {   cout << "library " << l << ": found " << dists[l].size( ) 
                        << " links, power = " << power[l] << "\n";    }
               if ( dists[l].size( ) >= 2 ) have_data = True;    }
          if ( !have_data )
          {    
               #pragma omp critical
               {    cout << "failed to find links for gap from " << m1
                         << " to " << m2 << ", old " << S.Gap(jg) << " +/- " 
                         << S.Dev(jg) << endl;
                    if (SHOW_START_STOP) cout << "stop " << m1 << endl;
                    fails++;    }
               continue;    }
          /*
          const double pow_mult = 3.0;
          for ( int l = 0; l < nlibs; l++ )
          {    if ( lfail[l] ) continue;
               if ( power[l] > pow_mult * Min(power) ) dists[l].clear( );    }
          */

          // Get start/stop points for last 10 kb of m1 and first 10 kb of m2.

          vec< vec<int> > starts1(nlibs), stops2(nlibs);
          for ( int l = 0; l < nlibs; l++ )
          {    if ( lfail[l] ) continue;
               for ( int jx = 0; jx < starts_fw[l][m1].isize( ); jx++ )
               {    int x = starts_fw[l][m1][jx];
                    if ( x <= tigs[m1].isize( ) - read_length[l]
                         && x >= tigs[m1].isize( ) - max_dist[l] && x >= 0 )
                    {    if ( uniq[m1][x] ) starts1[l].push_back(x);    }    }
               for ( int jx = 0; jx < stops_rc[l][m2].isize( ); jx++ )
               {    int x = stops_rc[l][m2][jx];
                    if ( x <= max_dist[l] - read_length[l] ) 
                         stops2[l].push_back(x);    }    }

          // Create a density map for the stops, smoothing.

          vec< vec<double> > density2(nlibs);
          for ( int l = 0; l < nlibs; l++ )
          {    if ( lfail[l] ) continue;
               density2[l].resize( max_dist[l], 0 );
               int smooth_delta = 50;
               for ( int j = 0; j < stops2[l].isize( ); j++ )
               {    for ( int x = stops2[l][j] - smooth_delta; 
                         x <= stops2[l][j] + smooth_delta; x++ )
                    {    if ( x >= 0 && x < max_dist[l] ) 
                              density2[l][x]++;    }    }    }
          for ( int l = 0; l < nlibs; l++ )
          {    if ( lfail[l] ) continue;
               double dm = Mean( density2[l] );
               // if ( tigs[m2].isize( ) < max_dist[l] )
               //      dm *= double(max_dist[l]) / double( tigs[m2].size( ) );
               for ( int j = 0; j < max_dist[l]; j++ )
                    density2[l][j] /= (dm*5.0);    }
     
          // Build (start1+x,start1,x) data.

          const int ndata = 1000000;
          vec< vec< triple<int,int,int> > > grunkle(nlibs);
          int rptr = 0;
          int r = rand[rptr];
          if ( ++rptr == rand.isize( ) ) rptr = 0;
          int ptr = r % D2vals[0].isize( ); // assume all D2vals of same size
          for ( int l = 0; l < nlibs; l++ )
          {    if ( lfail[l] ) continue;
               grunkle[l].resize(ndata);
               if ( starts1[l].empty( ) ) continue;
               for ( int j = 0; j < ndata; j++ )
               {    int r = rand[rptr];
                    if ( ++rptr == rand.isize( ) ) rptr = 0;
                    int start1 = starts1[l][ r % starts1[l].isize( ) ];
                    int x = D2vals[l][ptr];
                    if ( ++ptr == D2vals[l].isize( ) ) ptr = 0;
                    grunkle[l][j] = make_triple( start1+x, start1, x );    }
               Sort( grunkle[l] );    }

          // Report result for true gap.

          if ( VERBOSITY >= 2 ) PRINT4( m1, m2, tigs[m1].size( ), tigs[m2].size( ) );
          if ( VERBOSITY >= 2 && VALIDATE )
          {    double KS = 0;
               vec< vec< pair<int,IntDistribution> > > computed(nlibs);
               Bool OK = False;
               for ( int l = 0; l < nlibs; l++ )
               {    if ( dists[l].empty( ) ) continue;
                    if ( lfail[l] ) continue;
                    double KSl = GapToKS( computed[l], grunkle[l], tgap, rand, rptr, 
                         tigs[m1].size( ), tigs[m2].size( ), read_length[l], 
                         max_dist[l], K, uniq, m1, m2, density2[l], dists[l],
                         starts_fw[l][m1], stops_rc[l][m2], half_cov[l], DISTS[l] );
                    if ( KSl < 2 ) 
                    {    KS += KSl;
                         OK = True;    }    }
               if ( !OK ) KS = 2;
               cout << "\nKS for true gap " << tgap << " is " << KS << "\n";

               // Predict number of links.

               vec<double> exp_links( nlibs, 0 );
               for ( int l = 0; l < nlibs; l++ )
               {    for ( int j = 0; j < starts_fw[l][m1].isize( ); j++ )
                    {    int start = starts_fw[l][m1][j];
                         int d = DISTS[l][ randomx( ) % DISTS[l].size( ) ];
                         int stop = start + d - tigs[m1].isize( ) - tgap;
                         if ( stop < read_length[l] ) continue;
                         if ( stop >= tigs[m2].isize( ) ) continue;

                         // Compare to expected coverage on m2.

                         int delta = 100;
                         int count2 = 0;
                         for ( int k = 0; k < stops_rc[l][m2].isize( ); k++ )
                         {    if ( Abs( stop - stops_rc[l][m2][k] ) <= delta )
                                   count2++;    }
                         double exp2 = half_cov[l] * double( 2*delta + 1 );
                         double add = Min( 1.0, double(count2) / exp2 );
                         exp_links[l] += add;    }
                    cout << "expected links for library " << l << " = " 
                         << exp_links[l] << ", found " << dists[l].size( )
                         << endl;    }    }

          // Try various gap sizes.  Do multiple passes using half the data.

          const int npasses = 50;
          vec< vec< pair<int,IntDistribution> > > computed(nlibs);
          vec< vec< vec<int> > > dists_half(npasses);
          for ( int pass = 0; pass < npasses; pass++ )
          {    dists_half[pass].resize(nlibs);
               for ( int l = 0; l < nlibs; l++ )
               {    if ( lfail[l] ) continue;
                    for ( int i = 0; i < dists[l].isize( )/2; i++ )
                    {    int r = rand[rptr];
                         if ( ++rptr == rand.isize( ) ) rptr = 0;
                         dists_half[pass][l].push_back( 
                              dists[l][ r % dists[l].isize( ) ] );    }
                    Sort( dists_half[pass][l] );    }    }

          // Do two outer passes.  The first pass uses the lower bound, whereas the
          // second does not.  The second determines the standard deviation.

          vec< vec<float> > best_gaps(2);
          if ( VERBOSITY >= 2 ) cout << Date( ) << ": starting main loop" << endl;
          Bool fail = False;
          for (int opass = 1; opass <= 2; opass++ )
          {    int succeeds = 0;
               for ( int pass = 0; pass < npasses; pass++ )
               {    
                    // Find the best gap.
     
                    Bool pfail = False;
                    if ( VERBOSITY >= 2 ) 
                    {    cout << "\nopass = " << opass << ", npass = " << pass+1 
                              << "\n";    }
                    int delta_init = 1000;
                    vec< pair<double,int> > KS_gap;
                    KS_gap.push( -1, 0 );
                    if ( delta_init <= max_overlap ) KS_gap.push( -1, -delta_init );
                    if ( 2*delta_init <= max_overlap ) 
                         KS_gap.push( -1, -2*delta_init );
                    /*
                    for ( int d = delta_init*2; d < Max(max_dist); d += delta_init )
                         KS_gap.push( -1, d );
                    */
                    while(1)
                    {    for ( int i = 0; i < KS_gap.isize( ); i++ )
                         {    if ( KS_gap[i].first < 0 )
                              {    KS_gap[i].first = 0;
                                   int gap = KS_gap[i].second;
                                   Bool OK = False;
                                   for ( int l = 0; l < nlibs; l++ )
                                   {    if ( dists[l].size( ) <= 1 ) continue;
                                        if ( lfail[l] ) continue;
                                        double KS = GapToKS( computed[l], 
                                             grunkle[l], gap, rand, rptr, 
                                             tigs[m1].size( ), tigs[m2].size( ), 
                                             read_length[l], max_dist[l], K, uniq, 
                                             m1, m2, density2[l], 
                                             dists_half[pass][l],
                                             starts_fw[l][m1], stops_rc[l][m2],
                                             half_cov[l], DISTS[l] );
                                        if ( KS < 2 ) 
                                        {    KS_gap[i].first += KS;
                                             OK = True;    }    }
                                   if ( !OK ) KS_gap[i].first = 2;
                                   if ( VERBOSITY >= 2 )
                                   {    double KS = KS_gap[i].first;
                                        PRINT2( gap, KS );    }    }    }
                         Bool OK = False;
                         for ( int i = 0; i < KS_gap.isize( ); i++ )
                              if ( KS_gap[i].first < 2 ) OK = True;
                         if ( !OK )
                         {    pfail = True;
                              break;    }
                         Sort(KS_gap);
                         vec<int> bests, tried, to_try;
                         for ( int l = 0; l < KS_gap.isize( ); l++ )
                         {    if ( KS_gap[l].first == KS_gap[0].first )
                                   bests.push_back( KS_gap[l].second );
                              tried.push_back( KS_gap[l].second );    }
                         Sort(tried);
                         int ntried = tried.size( );
                         for ( int m = 0; m < bests.isize( ); m++ )
                         {    int x = bests[m];
                              int p = BinPosition( tried, x );
                              int delta = ( p == 0 ? delta_init : (x-tried[p-1])/2 );
                              if ( delta > 0 )
                              {    if ( opass == 1 )
                                   {    while ( -(x-delta) > max_overlap ) 
                                             delta /= 2;    }
                                   if ( delta > 0 ) 
                                        to_try.push_back( x - delta );    }
                              delta = ( p == ntried - 1 ? 
                                   delta_init : (tried[p+1]-x)/2 );
                              if ( delta > 0 ) to_try.push_back( x + delta );    }
                         if ( to_try.empty( ) ) break;
                         for ( int m = 0; m < to_try.isize( ); m++ )
                              KS_gap.push( -1, to_try[m] );    }
                    if ( !pfail ) succeeds++;
                    int best_gap = KS_gap[0].second;
                    double best_KS = KS_gap[0].first;
                    if ( VERBOSITY >= 2 ) 
                    {    cout << "gap = " << best_gap << ", KS = " << best_KS 
                              << " ***\n";    }
                    best_gaps[opass-1].push_back(best_gap);    }
               if ( succeeds == 0 ) fail = True;    }
          if ( VERBOSITY >= 2 ) cout << Date( ) << ": done with main loop" << endl;
          if (fail)
          {
               #pragma omp critical
               {    fails++;
                    cout << "failed to estimate gap from " << m1
                         << " to " << m2 << ", old " << S.Gap(jg) << " +/- " 
                         << S.Dev(jg) << endl;
                    if (SHOW_START_STOP) cout << "stop " << m1 << endl;    }
               continue;    }

          // Compute mean and standard deviation.

          NormalDistribution D1 = SafeMeanStdev(best_gaps[0]);
          NormalDistribution D2 = SafeMeanStdev(best_gaps[1]);
          int gap = int( round( D1.mu_ ) );
          int dev = int( round( D2.sigma_ / sqrt(2.0) ) );

          // Get alternate estimate of standard deviation.

          /*
          vec<double> ldev(nlibs);
          for ( int l = 0; l < nlibs; l++ )
          {    int count = 0;
               for ( int j = 0; j < dists[l].isize( ); j++ )
                    if ( dists[l][j] + gap < max_dist[l] ) count++;
               ldev[l] = lib_dev[l] / sqrt( double(count) );    }
          double zz = 0;
          for ( int l = 0; l < nlibs; l++ )
               zz += 1.0 / ( ldev[l] * ldev[l] );
          int alt_dev = int( ceil( 1.0 / sqrt(zz) ) );
          dev = Max( dev, alt_dev );
          */

          // Announce what's happening.

          if (SHOW_START_STOP)
          {
               #pragma omp critical
               {    cout << "stop " << m1 << endl;    }    }
          if ( VERBOSITY >= 2 ) cout << "\n";
          #pragma omp critical
          {    if ( VERBOSITY >= 1 )
               {    cout << m1 << "[" << tigs[m1].size( )  << "b] --> " << m2 << "[" 
                         << tigs[m2].size( ) << "b]: " << gap << " +/- " << dev;    }
               if (VALIDATE)
               {    double offby = double( Abs( gap - tgap ) ) / double(dev);
                    double aoffby 
                         = double( Abs( S.Gap(jg) - tgap ) ) / double( S.Dev(jg) );
                    if ( VERBOSITY >= 1 )
                    {    cout << ", true " << tgap << ", off " 
                              << setiosflags(ios::fixed) << setprecision(1)
                              << offby << " devs" << resetiosflags(ios::fixed);    }
                    results.push( gap, tgap, S.Gap(jg) );    
                    devs.push_back(dev);
                    adevs.push_back( S.Dev(jg) );
                    offbys.push_back(offby), aoffbys.push_back(aoffby);    }
               if ( VERBOSITY >= 1 )
               {    cout << ", old " << S.Gap(jg) << " +/- " << S.Dev(jg)
                         << endl;    }    }    }

     // Summarize results.

     if ( VALIDATE && ( TIGS == "" || tigs_to_process.size( ) > 1 ) )
     {    
          // Report count and devs.

          Sort(devs), Sort(adevs), Sort(offbys), Sort(aoffbys);
          int N = 25;
          int unassayed = 0;
          if ( tigs_to_process.empty( ) )
               unassayed = true_gap_computed.isize( ) - Sum(true_gap_computed);
          else
          {    for ( int j = 0; j < tigs_to_process.isize( ); j++ )
               {    if ( !true_gap_computed[ tigs_to_process[j] ] )
                         unassayed++;    }    }
          if ( devs.empty( ) )
          {    cout << "\ngaps:\n" << devs.size( ) << " + " << fails << " fails"
                    << "\n+ " << unassayed << " unassayed\n";
               cout << "No gaps computed, report ends here." << endl;
               return 0;    }
          ostringstream outw1, outw2;
          outw1 << "\ngaps:\n" << devs.size( ) << " + " << fails << " fails"
               << "\n+ " << unassayed << " unassayed\n";
          outw2 << "\ndevs:\n" << "new median dev = " << Median(devs) << "\n"
               << "old median dev = " << Median(adevs) << "\n";
          PrintSideBySide( outw1.str( ), outw2.str( ), N );
          double adj = double(Median(devs)) / double(Median(adevs));

          // Report off by in devs.  The last column old! is an adjusted version,
          // taking into account the difference between devs and adevs.

          cout << "\n";
          ostringstream outz1;
          outz1 << "percent off by:\n";
          vec< vec<String> > rows;
          vec<String> row;
          row.push_back( "dev", "new", "old", "old!" );
          rows.push_back(row);
          row.clear( );
          for ( double x = 3.0; x <= 5.0; x++ )
          {    int count = 0, acount = 0, aacount = 0;
               for ( int j = 0; j < offbys.isize( ); j++ )
               {    if ( offbys[j] > x ) count++;
                    if ( aoffbys[j] > x ) acount++;
                    if ( aoffbys[j] > x * adj ) aacount++;    }
               ostringstream out1, out2, out3;
               out1 << setiosflags(ios::fixed) << setprecision(1)
                    << 100.0 * double(count) / double( offbys.size( ) );
               out2 << setiosflags(ios::fixed) << setprecision(1)
                    << 100.0 * double(acount) / double( aoffbys.size( ) );
               out3 << setiosflags(ios::fixed) << setprecision(1)
                    << 100.0 * double(aacount) / double( aoffbys.size( ) );
               row.push_back( ToString(x), out1.str( ), out2.str( ), out3.str( ) );  
               rows.push_back(row);
               row.clear( );    }
          PrintTabular( outz1, rows, 2, "rrrr" );
          String z1 = outz1.str( );

          // Report absolute off by.

          int err_lt1 = 0;
          int err_lt10 = 0;
          int err_lt100 = 0;
          int err_lt1000 = 0;
          int err_lt10000 = 0;
          int err_huge = 0;
          int xerr_lt1 = 0;
          int xerr_lt10 = 0;
          int xerr_lt100 = 0;
          int xerr_lt1000 = 0;
          int xerr_lt10000 = 0;
          int xerr_huge = 0;
          for ( int j = 0; j < results.isize( ); j++ )
          {    int err = Abs( results[j].first - results[j].second );
               if ( err < 1 ) err_lt1++;
               else if ( err < 10 ) err_lt10++;
               else if ( err < 100 ) err_lt100++;
               else if ( err < 1000 ) err_lt1000++;
               else if ( err < 10000 ) err_lt10000++;
               else err_huge++;
               int xerr = Abs( results[j].third - results[j].second );
               if ( xerr < 1 ) xerr_lt1++;
               else if ( xerr < 10 ) xerr_lt10++;
               else if ( xerr < 100 ) xerr_lt100++;
               else if ( xerr < 1000 ) xerr_lt1000++;
               else if ( xerr < 10000 ) xerr_lt10000++;
               else xerr_huge++;    }
          rows.clear( );
          row.push_back( "error", "new", "old" );
          rows.push_back(row), row.clear( );
          row.push_back( "< 1", ToString(err_lt1), ToString(xerr_lt1) );
          rows.push_back(row), row.clear( );
          row.push_back( "< 10", ToString(err_lt10), ToString(xerr_lt10) );
          rows.push_back(row), row.clear( );
          row.push_back( "< 100", ToString(err_lt100), ToString(xerr_lt100) );
          rows.push_back(row), row.clear( );
          row.push_back( "< 1000", ToString(err_lt1000), ToString(xerr_lt1000) );
          rows.push_back(row), row.clear( );
          row.push_back( "< 10000", ToString(err_lt10000), ToString(xerr_lt10000) );
          rows.push_back(row), row.clear( );
          row.push_back( ">= 10000", ToString(err_huge), ToString(xerr_huge) );
          rows.push_back(row), row.clear( );
          ostringstream outz2;
          outz2 << "absolute error:\n";
          PrintTabular( outz2, rows, 2, "lrr" );
          String z2 = outz2.str( );

          // Merge tables.

          PrintSideBySide( z1, z2, N );    }

     // Report time used.

     cout << "\ntime used = " << TimeSince(clock) << " (" << getHostName()
             << ")\n\n";
     cout << "command: " << command.TheCommand( ) << "\n\n";    }
