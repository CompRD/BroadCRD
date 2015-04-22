///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// LongReadConsensus.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "PrintAlignment.h"
#include "Map.h"
#include "Set.h"
#include "graph/Digraph.h"
#include "math/HoInterval.h"
#include "paths/BigMapTools.h"
#include "paths/CorrectLongReadsTools2.h"
#include "paths/LongReadTools.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/UnibaseUtils.h"
#include "paths/Unipath.h"
#include "paths/Useq.h"
#include "paths/long/CreateGenome.h"
#include "random/Bernoulli.h"
#include "random/Random.h"

void Assess( const vec< vec<int> >& reads, const vecbasevector& genome,
     const vecbasevector& genome2, const vecbasevector& unibases, const int K, 
     const int LG, const VecIntPairVec& Glocs,
     const vec<int>& uperfect )
{    cout << "\n" << Date( ) << ": finding perfect reads" << endl;
     vec<Bool> lperfect( reads.size( ) );
     vec<placementy> places_all;
     vec< vec<ho_interval> > cov( genome.size( ) );
     int nz = 0;
     cout << "\n";
     for ( int id = 0; id < reads.isize( ); id++ )
     {    const vec<int>& x = reads[id];
          if ( x.empty( ) ) continue;
          nz++;
          basevector b = unibases[ x[0] ];
          for ( int j = 1; j < x.isize( ); j++ )
          {    b.resize( b.isize( ) - (K-1) );
               b = Cat( b, unibases[ x[j] ] );    }
          vec<placementy> places 
               = FindGenomicPlacementsY( 0, b, LG, genome2, Glocs );
          for ( int m = 0; m < places.isize( ); m++ )
          {    placementy p = places[m];
               cov[p.g].push( p.pos, p.Pos );
               int ng = genome2[p.g].isize( ) / 2;
               if ( p.pos >= ng ) continue;    }
          places_all.append(places);
          static int count(0);
          int pc = 200;
          if ( places.empty( ) && ++count % pc == 0 )
          {    cout << "[" << count/pc << "] " << id << " imperfect";
               for ( int j = 0; j < x.isize( ); j++ )
               {    cout << " " << x[j];
                    if ( uperfect[ x[j] ] == 0 ) cout << "[FALSE]";    }
               cout << "\n";    }
          lperfect[id] = places.nonempty( );    }
     cout << "\n";
     PRINT(nz);
     cout << PERCENT_RATIO( 3, Sum(lperfect), nz ) << " of reads are "
          << "perfect" << endl;
     vec< vec<ho_interval> > cov2( genome.size( ) ), cov3( genome.size( ) );
     for ( int g = 0; g < (int) genome.size( ); g++ )
     {    Sort( cov[g] );
          for ( int j = 0; j < cov[g].isize( ); j++ )
          {    int k;
               for ( k = j + 1; k < cov[g].isize( ); k++ )
                    if ( cov[g][k].Start( ) > cov[g][j].Start( ) ) break;
               int stop = 0;
               for ( int l = j; l < k; l++ )
                    stop = Max( stop, cov[g][l].Stop( ) );
               ho_interval h( cov[g][j].Start( ), stop );
               cov2[g].push_back(h);
               j = k - 1;    }    }
     for ( int g = 0; g < (int) genome.size( ); g++ )
     {    for ( int j = 0; j < cov2[g].isize( ); j++ )
          {    int k;
               for ( k = j + 1; k < cov2[g].isize( ); k++ )
                    if ( cov2[g][k].Stop( ) > cov2[g][j].Stop( ) ) break;
               cov3[g].push_back( cov2[g][j] );
               j = k - 1;    }    }
     vec<int> overlaps;
     Bool verbose = False;
     for ( int g = 0; g < (int) genome.size( ); g++ )
     {    if (verbose) cout << "\ncoverage of reference " << g << endl;
          for ( int j = 0; j < cov3[g].isize( ); j++ )
          {    if (verbose) cout << cov3[g][j];
               if ( j > 0 ) 
               {    int overlap = cov3[g][j-1].Stop( ) - cov3[g][j].Start( );
                    if (verbose) cout << " [" << overlap << "]";
                    if ( overlap < 1000 ) overlaps.push_back(overlap);    }
               if (verbose) cout << "\n";    
               if ( cov3[g][j].Stop( ) > genome[g].isize( ) + 1000 ) break;    }    }
     int lt0 = 0, lt200 = 0, lt500 = 0, lt640 = 0, lt1000 = 0;
     for ( int i = 0; i < overlaps.isize( ); i++ )
     {    if ( overlaps[i] < 0 ) lt0++;
          else if ( overlaps[i] < 200 ) lt200++;
          else if ( overlaps[i] < 500 ) lt500++;
          else if ( overlaps[i] < 640 ) lt640++;
          else lt1000++;    }
     PRINT5( lt0, lt200, lt500, lt640, lt1000 );    }

// Count the number of true locations that we should have found.
// Note that there's a bug - we don't correctly handle the case where
// the unibases wrap around the beginning of the genome.

int TrueLocs( const int test_id, const int left_flank, const int right_flank, 
     const int K, const vecbasevector& unibases, const vec<int>& nkmers,
     const vec< vec<int> >& nexts, const vec<int>& to_rc, 
     const vecbasevector& genome, const vecbasevector& genome2, const int LG,
     const VecIntPairVec& Glocs )
{
     int gcount = 0;
     vec<placementy> places = FindGenomicPlacementsY( 
          0, unibases[test_id], LG, genome2, Glocs );
     for ( int j = 0; j < places.isize( ); j++ )
     {    if ( places[j].pos >= genome[ places[j].g ].isize( ) ) continue;
          int left_ext = 0, right_ext = 0;
          int u = test_id;
          placementy p_curr = places[j];
          while( right_ext < right_flank )
          {    for ( int l = 0; l < nexts[u].isize( ); l++ )
               {    int v = nexts[u][l];
                    vec<placementy> places_v = FindGenomicPlacementsY( 
                         0, unibases[v], LG, genome2, Glocs );
                    for ( int r = 0; r < places_v.isize( ); r++ )
                    {    placementy p_next = places_v[r];
                         if ( p_next.g != p_curr.g ) continue;
                         if ( p_next.fw != p_curr.fw ) continue;
                         if ( ( p_curr.fw && p_curr.Pos == p_next.pos + K - 1 )
                              || ( !p_curr.fw && p_next.Pos == p_curr.pos + K - 1 ) )
                         {    right_ext += nkmers[v];
                              u = v;
                              p_curr = p_next;
                              goto next_right;    }    }    }
               break;
               next_right: continue;    }
          if ( right_ext < right_flank ) continue;
          u = test_id;
          p_curr = places[j];
          while( left_ext < left_flank )
          {    for ( int l = 0; l < nexts[ to_rc[u] ].isize( ); l++ )
               {    int v = to_rc[ nexts[ to_rc[u] ][l] ];
                    vec<placementy> places_v = FindGenomicPlacementsY( 
                         0, unibases[v], LG, genome2, Glocs );
                    for ( int r = 0; r < places_v.isize( ); r++ )
                    {    placementy p_next = places_v[r];
                         if ( p_next.g != p_curr.g ) continue;
                         if ( p_next.fw != p_curr.fw ) continue;
                         if ( ( p_curr.fw && p_next.Pos == p_curr.pos + K - 1 )
                              || ( !p_curr.fw && p_curr.Pos == p_next.pos + K - 1 ) )
                         {    left_ext += nkmers[v];
                              u = v;
                              p_curr = p_next;
                              goto next_left;    }    }    }
               break;
               next_left: continue;    }
          if ( left_ext >= left_flank ) gcount++;    }
     return gcount;    }

Bool TrimByFilledFragments( vec<int>& x, const vec< vec<int> >& nexts,
     const vec< vec<int> >& nexts_count,
     vec< StdMap< vec<int>, vec< triple< int, int, int > > > >& segmap,
     const Bool trim )
{
     // Define heuristics.

     const double p_max = 0.05;
     const int k = 10;
     const int max_score = 4;
     double p = 1.0/double(k+1);

     // Trim ucores containing branches that are weak according to the filled
     // fragments.  

     int left_trim = 0, right_trim = 0;
     for ( int j = 0; j < x.isize( ) - 1; j++ )
     {    int u = x[j], v = x[j+1], score = 0; 
          int best_score = Sum( nexts_count[u] );
          for ( int l = 0; l < nexts[u].isize( ); l++ )
               if ( nexts[u][l] == v ) score = nexts_count[u][l];

          // Statistical test, see below.

          int n = best_score, r = score;
          if ( score <= max_score && n >= 1 && BinomialSum(n, r, p) <= p_max ) 
          {    int left = j + 1, right = x.isize( ) - (j+2) + 1;
               if ( left <= left_trim || right <= right_trim ) continue;
               if ( left - left_trim <= right - right_trim ) left_trim = left;
               else right_trim = right;    }    }
     if ( !trim && ( left_trim > 0 || right_trim > 0 ) ) return False;
     x.SetToSubOf( x, left_trim, x.isize( ) - left_trim - right_trim );

     // Find segments having size between 3 and 6 that appear with low frequency in 
     // the filled fragments, and kill the ucores that contain them.  For example, 
     // for size 4, wherever we find a b c d in the ucores, we let 
     // best_score = max mult(a',b,c,d') in the filled fragments (a', b' varying)
     //      score =     mult(a, b,c,d ) in the filled fragments,
     // and if best_score >= 10 and best_score > 10 * score, we trim the ucore.

     int max_seg = segmap.isize( ) - 1;
     for ( int segsize = 3; segsize <= max_seg; segsize++ )
     {    int left_trim = 0, right_trim = 0;
          for ( int j = 0; j <= x.isize( ) - segsize; j++ )
          {    vec<int> key;
               key.SetToSubOf( x, j+1, segsize-2 );

               // add this if you want to parallelize the loop
               // if ( segmap[segsize].find(key) == segmap[segsize].end( ) ) 
               //      continue;

               vec< triple<int,int,int> > value = segmap[segsize][key];

               int best_score = 0, score = 0;
               for ( int l = 0; l < value.isize( ); l++ )
               {    best_score += value[l].third;
                    // best_score = Max( best_score, value[l].third );
                    if ( value[l].first == x[j] 
                         && value[l].second == x[ j + segsize - 1 ] )
                    {    score = value[l].third;    }    }

               // Suppose that the actual prevalence of the given site is 
               // k-fold lower than the most prevalent site (taking k = 5).
               // Then we require that the probability of seeing r or less
               // observations (out of n) is at most p_max = 5%.  Thus we test
               // for sum_{y=0)^r (n choose y) p^y (1-p)^(n-y) <= p_max.

               int n = best_score, r = score;
               if ( score <= max_score && n >= 1 && BinomialSum(n, r, p) <= p_max ) 
               {    int left = j + 1, right = x.isize( ) - (j+segsize) + 1;
                    if ( left <= left_trim || right <= right_trim ) continue;
                    if ( left - left_trim <= right - right_trim ) left_trim = left;
                    else right_trim = right;    }    }
          if ( !trim && ( left_trim > 0 || right_trim > 0 ) ) return False;
          x.SetToSubOf(x, left_trim, x.isize( ) - left_trim - right_trim);    }
     return True;    }

Bool Match( const vec<int>& x1, const int p1, const vec<int>& x2, const int p2 )
{    for ( int j1 = p1; j1 >= 0; j1-- )
     {    int j2 = j1 + p2 - p1;
          if ( j2 < 0 ) break;
          if ( x1[j1] != x2[j2] ) return False;    }
     for ( int j1 = p1 + 1; j1 < x1.isize( ); j1++ )
     {    int j2 = j1 + p2 - p1;
          if ( j2 == x2.isize( ) ) break;
          if ( x1[j1] != x2[j2] ) return False;    }
     return True;    }

// Coverage: compute the long read coverage of a chain of unibases, each overlapping
// the next by K-1.  As computed here, for a read to cover at all, it must cover the
// entire chain of unipaths, including the ends of the ends.

void Coverage( vec<int> u, const int K, const int L, const vecbasevector& unibases, 
     const vec<int>& to_rc, const vec< vec<int> >& cores, 
     const vec< vec<int> >& ucores, const vec< vec< pair<int,int> > >& cindex,
     const vec< vec< pair< int, vec< pair<int,int> > > > >& aligns, vec<int>& cov )
{
     int nbases = unibases[ u[0] ].size( );
     for ( int j = 1; j < u.isize( ); j++ )
          nbases += unibases[ u[j] ].isize( ) - (K-1);
     cov.resize( nbases - (L-1), 0 );
     for ( int pass = 1; pass <= 2; pass++ )
     {    if ( pass == 2 )
          {    u.ReverseMe( );
               for ( int j = 0; j < u.isize( ); j++ )
                    u[j] = to_rc[ u[j] ];    }
          for ( int m = 0; m < cindex[ u[0] ].isize( ); m++ )
          {    int id = cindex[ u[0] ][m].first, p = cindex[ u[0] ][m].second;
               if ( id >= cores.isize( ) / 2 ) continue;
               const vec<int> &x = ucores[id], &y = cores[id];
               // if ( p + u.isize( ) > x.isize( ) ) continue;
               if ( p == 0 || p + u.isize( ) > x.isize( ) - 1 ) continue;
               Bool mismatch = False;
               for ( int j = 1; j < u.isize( ); j++ )
               {    if ( u[j] != x[p+j] )
                    {    mismatch = True;
                         break;    }    }
               if (mismatch) continue;
               int pos = 0;
               for ( int j = 0; j < u.isize( ); j++ )
               {    int n = y[p+j];
                    const vec< pair<int,int> >& P = aligns[id][n].second;
                    for ( int z = 0; z < P.isize( ); z++ )
                    {    int upos = P[z].second;
                         if ( j < u.isize( ) - 1 
                              && upos >= unibases[ u[j] ].isize( ) - (K-1) )
                         {    continue;    }
                         int cp = pos + upos;
                         if ( pass == 2 ) cp = nbases - cp - L;
                         cov[cp]++;    }
                    pos += unibases[ u[j] ].isize( ) - (K-1);    }    }    }    }

void FindIllegalPaths( const int len_max, const int K, const int L, 
     const vecbasevector& unibases, const vec<int>& nkmers, const vec<int>& to_rc, 
     const vec< vec<int> >& cores, const vec< vec<int> >& ucores, 
     const vec< vec< pair< int, vec< pair<int,int> > > > >& aligns, 
     const vecbasevector& genome2, const int LG, 
     const VecIntPairVec& Glocs, vec< vec< vec<int> > >& illegal,
     const Bool VERBOSE )
{
     // Compute cindex.

     cout << Date( ) << ": compute cindex" << endl;
     int nuni = unibases.size( );
     vec< vec< pair<int,int> > > cindex(nuni);
     for ( int id = 0; id < ucores.isize( ); id++ )
     for ( int p = 0; p < ucores[id].isize( ); p++ )
          cindex[ ucores[id][p] ].push( id, p );

     // Get read length distribution, computed in a certain way.

     vec<int> lens;
     for ( int id = 0; id < ucores.isize( )/2; id++ )
     {    const vec<int>& x = ucores[id];
          int nk = 0;
          for ( int j = 1; j < x.isize( ) - 1; j++ )
               nk += nkmers[ x[j] ];
          if ( nk > 0 ) lens.push_back(nk);    }
     Sort(lens);

     // Compute expected coverage.

     const int len_low = 300;
     const int len_high = 500;
     vec<double> cove;
     for ( int u = 0; u < nuni; u++ )
     {    if ( nkmers[u] < len_low || nkmers[u] > len_high ) continue;
          vec<int> x;
          x.push_back(u);
          vec<int> cov;
          Coverage( x, K, L, unibases, to_rc, cores, ucores, cindex, aligns, cov );
          double mean_cov = double( Sum(cov) ) / double( cov.size( ) );
          cove.push_back(mean_cov);    }
     Sort(cove);
     // Note danger if not enough values *******************************************
     double mcov = cove.size() > 0 ? cove[ cove.size( ) / 2 ] : 0;
     if ( cove.size() == 0 )
       cout << "\n        WARNING!!!!!!!!!!WARNING!!!!!!!!!!WARNING!!!!!!!!\n  found no long read coverage on unibases" << endl;

     // Look for improbable unipath sequences.

     const double cov_mult = 0.4;
     const int max_len = 1000;
     vec<int> rands;
     const int rsize = 1000000;
     for ( int i = 0; i < rsize; i++ )
          rands.push_back( randomx( ) );
     vec< StdSet< vec<int> > > X( len_max + 1 );
     illegal.resize( len_max + 1 );
     for ( int len = 1; len <= len_max; len++ )
     {    
          // Get candidates.

          vec< vec<int> > candidates;
          for ( int id = 0; id < ucores.isize( ); id++ )
          for ( int p = 0; p <= ucores[id].isize( ) - len; p++ )
          {    vec<int> x;
               x.SetToSubOf( ucores[id], p, len );
               if ( Member( X[len], x ) ) continue;
               X[len].insert(x);
               Bool ill = False;
               for ( int l = 1; l < len; l++ )
               {    if (ill) break;
                    for ( int j = 0; j < illegal[l].isize( ); j++ )
                    {    if (ill) break;
                         if ( x.Contains( illegal[l][j] ) )
                              ill = True;    }    }
               if (ill) continue;
               candidates.push_back(x);    }

          // Process candidates.

          #pragma omp parallel for
          for ( int ci = 0; ci < candidates.isize( ); ci++ )
          {    const vec<int>& x = candidates[ci];
               int nk = 0;
               for ( int j = 0; j < x.isize( ); j++ )
                    nk += nkmers[ x[j] ];
               int lenf = LowerBound( lens, nk );
               double len_frac = double(lenf) / double( lens.size( ) );
               if ( nk > max_len ) continue;
               vec<int> cov;
               Coverage( x, K, L, unibases, to_rc, cores, ucores, 
                    cindex, aligns, cov );
               int min_cov = Min(cov);
               double mean_cov = double( Sum(cov) ) / double( cov.size( ) );
               int z = 0;
               for ( int j = 0; j < cov.isize( ); j++ )
                    if ( cov[j] == 0 ) z++;

               // Initial test to avoid wasting time.

               double cov_bound = cov_mult * mcov * (1.0-len_frac);
               if ( mean_cov >= cov_bound ) continue;

               // Estimate the probability that a random read with 25% error rate
               // would fail to align to this.

               int rand_ptr = 0, n = nk + K - 1;
               vec<ho_interval> required;
               if ( x.solo( ) ) required.push( 0, n );
               else
               {    int pos = 0;
                    for ( int j = 0; j < x.isize( ) - 1; j++ )
                    {    pos += nkmers[ x[j] ];
                         required.push( pos, pos + K - 1 );    }    }
               const int sim_count = 200;
               int fails = 0;
               for ( int j = 0; j < sim_count; j++ )
               {    vec<int> mis( n/4 );
                    for ( int r = 0; r < n/4; r++ )
                    {    int rx = rands[rand_ptr++];
                         if ( rand_ptr == rsize ) rand_ptr = 0;
                         mis[r] = rx % n;    }
                    UniqueSort(mis);
                    vec<ho_interval> goods;
                    if ( mis.empty( ) ) goods.push( 0, n );
                    else
                    {    if ( mis[0] >= L ) goods.push( 0, mis[0] );
                         for ( int j = 1; j < mis.isize( ); j++ )
                         {    if ( mis[j] - mis[j-1] >= L )
                                   goods.push( mis[j-1], mis[j] );    }
                         if ( n - mis.back( ) >= L )
                              goods.push( mis.back( ), n );    }
                    for ( int r = 0; r < required.isize( ); r++ )
                    {    const ho_interval& h = required[r];
                         Bool ok = False;
                         for ( int i = 0; i < goods.isize( ); i++ )
                         {    if ( Overlap( h, goods[i] ) >= L )
                              {    ok = True;
                                   break;    }    }
                         if ( !ok )
                         {    fails++;
                              break;    }    }    }

               double fail_frac = double(fails) / double(sim_count);
               mean_cov *= 1.0 / ( 1.0 - fail_frac );
               cov_bound = cov_mult * mcov * (1.0-len_frac);

               if ( mean_cov < cov_bound )
               {    
                    #pragma omp critical
                    {    illegal[len].push_back(x);    }
                    ostringstream hout;
                    hout << "\ncoverage of unibase seq";
                    for ( int j = 0; j < x.isize( ); j++ )
                         hout << " " << x[j];
                    hout << " (";
                    for ( int j = 0; j < x.isize( ); j++ )
                    {    if ( j > 0 ) hout << "+";
                         hout << nkmers[ x[j] ];    }
                    hout << "=" << nk << " kmers, " << x.size( ) << " units)\n";
                    for ( int j = 0; j < cov.isize( ); j++ )
                    {    if ( j > 0 ) hout << " ";
                         hout << cov[j];    }
                    hout << "\n";
                    PRINT2_TO( hout, mean_cov, cov_bound );
                    PRINT2_TO( hout, z, mean_cov );    
                    basevector u = unibases[ x[0] ];
                    for ( int j = 1; j < x.isize( ); j++ )
                    {    u.resize( u.isize( ) - (K-1) );
                         u = Cat( u, unibases[ x[j] ] );    }
                    if ( genome2.size( ) > 0 )
                    {    vec<placementy> places = FindGenomicPlacementsY( 
                              0, u, LG, genome2, Glocs );
                         if ( places.nonempty( ) ) 
                              hout << "WARNING: VALID!" << endl;    }
                    if (VERBOSE)
                    {    
                         #pragma omp critical
                         {    cout << hout.str( );    }    }    }    }    }    }

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_Int(K);
     CommandArgument_Int(KBIG);
     CommandArgument_Bool_OrDefault(WRITE, False);
     CommandArgument_Int_OrDefault_Doc(TRACE, -1, "search for this unibase");
     CommandArgument_Int_OrDefault_Doc(TEST_ID, -1, "only extend this unibase");
     CommandArgument_Int_OrDefault(MC_REL, 20);
     CommandArgument_Int_OrDefault(MC_BOUND, 20);
     CommandArgument_Int_OrDefault(MAX_RADIUS, 600);
     CommandArgument_Int_OrDefault(MIN_SEED, 30);
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0,
       "Number of threads to use (use all available processors if set to 0)");
     CommandArgument_Bool_OrDefault(VALIDATE, False);
     CommandArgument_Bool_OrDefault(VERBOSE, False);
     CommandArgument_Bool_OrDefault(EXPERIMENT, False);
     CommandArgument_Bool_OrDefault(EXPERIMENT_EDIT, False);
     CommandArgument_Int_OrDefault(EXPERIMENT_TEST, -1);
     CommandArgument_Int_OrDefault(EXPERIMENT_VERBOSITY, 0);
     CommandArgument_Bool_OrDefault(EXP2, True);
     CommandArgument_Bool_OrDefault(EXP2_SAVE, True);
     EndCommandArguments;

     // Define directories.

     cout << Date( ) << ": begin" << endl;
     String data_dir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String IN_HEAD = "extended";

     // Thread control.
     
     NUM_THREADS = configNumThreads(NUM_THREADS); 
     omp_set_num_threads( NUM_THREADS );

     // Load and hash genome.

     vecbasevector genome, genome2;
     const int LG = 12;
     VecIntPairVec Glocs;
     int ng;
     if (VALIDATE)
     {    genome.ReadAll( data_dir + "/genome.fastb" );
          ng = genome.size( );
          genome2.resize( genome.size( ) ); 
          for ( size_t j = 0; j < genome.size( ); j++ )
               genome2[j] = Cat( genome[j], genome[j] );
          CreateGlocs( genome2, LG, Glocs );    }

     // Load unibases and compute ancillary data.

     cout << Date( ) << ": load unibases" << endl;
     vecbasevector unibases( run_dir + "/" + IN_HEAD + ".unibases.k" + ToString(K) );
     int nuni = unibases.size( );
     PRINT(nuni);
     vec<int> to_rc, nkmers(nuni);
     UnibaseInvolution( unibases, to_rc );
     vec< vec<int> > nexts;
     GetNexts( K, unibases, nexts );
     for ( int u = 0; u < nuni; u++ )
          nkmers[u] = unibases[u].isize( ) - (K-1);

     // Determine which unibases are perfect.

     vec<int> uperfect( nuni, 0 );
     int nuperfect = 0;
     if (VALIDATE)
     {    cout << Date( ) << ": finding perfect unibases" << endl;
          for ( int u = 0; u < nuni; u++ )
          {    vec<placementy> places
                    = FindGenomicPlacementsY( 0, unibases[u], LG, genome2, Glocs );
               for ( int j = 0; j < places.isize( ); j++ )
               {    if ( places[j].pos < genome[ places[j].g ].isize( ) ) 
                         uperfect[u]++;    }
               if ( uperfect[u] > 0 ) nuperfect++;    }
          cout << PERCENT_RATIO( 3, nuperfect, nuni ) << " of unibases are "
               << "perfect" << endl;    }

     // Load filled fragments and map them back to the unibases.
     
     vec< vec<int> > fillseqs;
     {
     ForceAssertEq( K, 96 );
     const int K = 96;
     vecbasevector fills( run_dir + "/filled_reads.fastb" );
     vec< triple<kmer<K>,int,int> > kmers_plus;
     vec<int64_t> starts;
     starts.push_back(0);
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    const basevector& u = unibases[i];
          starts.push_back( starts.back( ) + Max( 0, u.isize( ) - K + 1 ) );    }
     kmers_plus.resize( starts.back( ) );
     #pragma omp parallel for
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    const basevector& u = unibases[i];
          kmer<K> x;
          for ( int j = 0; j <= u.isize( ) - K; j++ )
          {    int64_t r = starts[i] + j;
               x.SetToSubOf( u, j );
               kmers_plus[r].first = x;
               kmers_plus[r].second = i;
               kmers_plus[r].third = j;    }    }
     ParallelSort(kmers_plus);
     vec< kmer<K> > kmers( kmers_plus.size( ) );
     for ( size_t i = 0; i < kmers.size( ); i++ )
          kmers[i] = kmers_plus[i].first;
     vec< vec<int> > fillx( fills.size( ) );
     PRINT( fills.size( ) );
     #pragma omp parallel for
     for ( size_t id = 0; id < fills.size( ); id++ )
     {    const basevector& f = fills[id];
          if ( f.isize( ) <= K ) continue;
          kmer<K> x;
          x.SetToSubOf( f, 0 );
          int64_t p = BinPosition( kmers, x );
          if ( p < 0 ) continue;
          int u = kmers_plus[p].second, pos = kmers_plus[p].third;
          if ( pos + f.isize( ) <= unibases[u].isize( ) ) continue;
          vec< pair<int,int> > locs;
          locs.push( u, pos );
          for ( int j = 1; j <= f.isize( ) - K; j++ )
          {    x.SetToSubOf( f, j );
               p = BinPosition( kmers, x );
               if ( p < 0 ) break;
               u = kmers_plus[p].second, pos = kmers_plus[p].third;
               locs.push( u, pos );    }
          if ( locs.isize( ) < f.isize( ) - K + 1 ) continue;
          vec<int> us;
          Bool fail = False;
          for ( int i = 0; i < locs.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < locs.isize( ); j++ )
                    if ( locs[j].first != locs[i].first ) break;
               if ( i > 0 && locs[i].second != 0 ) fail = True;
               if ( i > 0 && locs[i-1].second
                    != unibases[ locs[i-1].first ].isize( ) - K )
               {    fail = True;    }
               for ( int k = i+1; k < j; k++ )
                    if ( locs[k].second != locs[k-1].second + 1 ) fail = True;
               if (fail) break;
               us.push_back( locs[i].first );
               i = j - 1;    }
          if (fail) continue;
          #pragma omp critical
          {    fillseqs.push_back(us);    }    }
     int nf = fillseqs.size( );
     for ( int i = 0; i < nf; i++ )
     {    vec<int> x = fillseqs[i];
          x.ReverseMe( );
          for ( int j = 0; j < x.isize( ); j++ )
               x[j] = to_rc[ x[j] ];
          fillseqs.push_back(x);    }
     PRINT( fillseqs.size( ) );
     Bool fverbose = False;
     if (fverbose)
     {    for ( int j = 0; j < fillseqs.isize( ); j++ )
          {    cout << "[" << j << " fill]";
               const vec<int>& x = fillseqs[j];
               for ( int i = 0; i < x.isize( ); i++ )
                    cout << " " << x[i];
               cout << "\n";    }    }

     if ( TRACE >= 0 )
     {    cout << "\nTRACE, filled fragments\n" << endl;
          for ( int j = 0; j < fillseqs.isize( ); j++ ) 
          {    const vec<int>& x = fillseqs[j]; 
               if ( !Member( x, TRACE ) ) continue; 
               cout << "[" << j << " fill]"; 
               for ( int i = 0; i < x.isize( ); i++ ) 
                    cout << " " << x[i]; 
               cout << "\n";    }     }
     }

     // Assess branches using filled fragments.

     cout << Date( ) << ": creating ffpairs" << endl;
     vec< pair<int,int> > ffpairs;
     for ( int i = 0; i < fillseqs.isize( ); i++ )
     {    const vec<int>& x = fillseqs[i];
          for ( int j = 0; j < x.isize( ) - 1; j++ )
               ffpairs.push( x[j], x[j+1] );    }
     cout << Date( ) << ": sorting" << endl;
     ParallelSort(ffpairs);
     cout << Date( ) << ": cataloging" << endl;
     vec< vec<int> > nexts_count(nuni);
     for ( int u = 0; u < nuni; u++ )
          nexts_count[u].resize( nexts[u].size( ), 0 );
     for ( int i = 0; i < ffpairs.isize( ); i++ )
     {    int j = ffpairs.NextDiff(i);
          int u = ffpairs[i].first, v = ffpairs[i].second;
          for ( int l = 0; l < nexts[u].isize( ); l++ )
               if ( nexts[u][l] == v ) nexts_count[u][l] = j - i;
          i = j - 1;    }

     // Load alignments.
     
     cout << Date( ) << ": loading alignments" << endl;
     vec< vec< pair< int, vec< pair<int,int> > > > > aligns;
     BinaryReader::readFile( run_dir + "/extended.long.read_aligns", &aligns );

     // Load paths.

     vec< vec< vec<int> > > paths;
     cout << Date( ) << ": loading paths" << endl;
     BinaryReader::readFile( run_dir + "/extended.long.read_paths", &paths );
     int nr = paths.size( );
     paths.append(paths);
     for ( int id = nr; id < 2*nr; id++ )
     {    for ( int m = 0; m < paths[id].isize( ); m++ )
               paths[id][m].ReverseMe( );    }
     
     // Load graphs.

     vec< digraphVE<int,int> > H;
     cout << Date( ) << ": loading H" << endl;
     BinaryReader::readFile( run_dir + "/extended.long.read_graph", &H );
     vec< digraphVE<int,int> > Hrc(H);
     for ( int i = 0; i < H.isize( ); i++ )
     {    Hrc[i].Reverse( );
          for ( int j = 0; j < Hrc[i].N( ); j++ )
               Hrc[i].VertMutable(j) = to_rc[ H[i].Vert(j) ];    }

     // Translate paths.

     vec< vec< vec<int> > > upaths(paths);
     for ( int id = 0; id < nr; id++ )
     {    for ( int m = 0; m < upaths[id].isize( ); m++ )
          {    vec<int>& x = upaths[id][m];
               for ( int j = 0; j < x.isize( ); j++ )
                    x[j] = H[id].Vert( x[j] );    }    }
     for ( int id = nr; id < 2*nr; id++ )
     {    for ( int m = 0; m < upaths[id].isize( ); m++ )
          {    vec<int>& x = upaths[id][m];
               for ( int j = 0; j < x.isize( ); j++ )
                    x[j] = Hrc[id-nr].Vert( x[j] );    }    }

     // Create index.

     vec< vec< triple<int,int,int> > > xindex(nuni);
     for ( int id = 0; id < 2*nr; id++ )
     {    for ( int m = 0; m < upaths[id].isize( ); m++ )
          {    const vec<int>& x = upaths[id][m];
               for ( int j = 0; j < x.isize( ); j++ )
                    xindex[ x[j] ].push( id, m, j );    }    }

     // Create cores.
     
     vec< vec<int> > cores(nr);
     for ( int id = 0; id < nr; id++ )
     {    if ( paths[id].empty( ) ) continue;
     
          // Find common substring.
          
          const vec<int>& x = paths[id][0];
          int left = -1, right = -1;
          Bool found = False;
          vec<int> common;
          for ( int trim = 0; trim <= x.isize( ); trim++ )
          {    if (found) break;
               for ( int j1 = 0; j1 <= trim; j1++ )
               {    if (found) break;
                    int j2 = trim - j1;
                    vec<int> y;
                    y.SetToSubOf( x, j1, x.isize( ) - trim );
                    Bool bad = False;
                    for ( int l = 0; l < paths[id].isize( ); l++ )
                    {    if ( !paths[id][l].Contains(y) )
                         {    bad = True;
                              break;    }    }
                    if ( !bad )
                    {    left = j1, right = j2;
                         common = y;
                         found = True;    }    }    }
          cores[id] = common;     }
     cores.append(cores);
     for ( int id = nr; id < 2*nr; id++ )
          cores[id].ReverseMe( );

     // Generate ucores.

     vec< vec<int> > ucores = cores;
     for ( int id = 0; id < nr; id++ )
     {    vec<int> x = cores[id];
          for ( int j = 0; j < x.isize( ); j++ )
               x[j] = H[id].Vert( x[j] );
          ucores[id] = x;    }
     for ( int id = nr; id < 2*nr; id++ )
     {    vec<int> x = cores[id];
          for ( int j = 0; j < x.isize( ); j++ )
               x[j] = Hrc[id-nr].Vert( x[j] );
          ucores[id] = x;    }
     vec< vec< pair<int,int> > > index(nuni);
     for ( int id = 0; id < ucores.isize( ); id++ )
     {    const vec<int>& x = ucores[id];
          for ( int j = 0; j < x.isize( ); j++ )
               index[ x[j] ].push( id, j );    }

     // Determine illegal sequences.

     vec< vec< vec<int> > > illegal;
     for ( int megapass = 1; megapass <= 2; megapass++ )
     {    cout << "\nstarting megapass " << megapass << endl;;
          int len_max = ( megapass == 1 ? 1 : 20 );
          const int L = 11;
          FindIllegalPaths( len_max, K, L, unibases, nkmers, to_rc, cores, ucores, 
               aligns, genome2, LG, Glocs, illegal, VERBOSE );

          // On the first pass we attempt to fix the cores.

          if ( megapass == 1 )
          {    vec<int> illegal1;
               for ( int j = 0; j < illegal[1].isize( ); j++ )
                    illegal1.push_back( illegal[1][j][0] );
               for ( int id = 0; id < ucores.isize( ); id++ )
               {    vec<int> &x = ucores[id], &y = cores[id];

                    // Prune until ends are legal.

                    while( x.nonempty( ) && Member( illegal1, x.front( ) ) )
                    {    x.erase( x.begin( ) );
                         y.erase( y.begin( ) );    }
                    while( x.nonempty( ) && Member( illegal1, x.back( ) ) )
                    {    x.resize( x.isize( ) - 1 );
                         y.resize( y.isize( ) - 1 );    }    }    }    }

     // Trace.

     if ( TRACE >= 0 )
     {    cout << "\nTRACE before cleaning bubbles\n" << endl;
          for ( int id = 0; id < ucores.isize( ); id++ ) 
          {    vec<int>& x = ucores[id]; 
               if ( !Member( x, TRACE ) ) continue;
               cout << "id = " << id << ", x =";
               for ( int j = 0; j < x.isize( ); j++ ) 
                         cout << " " << x[j]; 
               cout << endl;    }    } 

     // Clean bubbles.

     vec< pair< vec<int>, vec<int> > > edits;
     CleanBubbles( K, unibases, nexts, upaths, cores, ucores, index, genome2, 
          LG, Glocs, H, Hrc, uperfect, VERBOSE );

     // Trace.

     if ( TRACE >= 0 )
     {    cout << "\nTRACE before filled fragment trimming\n" << endl;
          for ( int id = 0; id < ucores.isize( ); id++ ) 
          {    vec<int>& x = ucores[id]; 
               if ( !Member( x, TRACE ) ) continue;
               cout << "id = " << id << ", x =";
               for ( int j = 0; j < x.isize( ); j++ ) 
                    cout << " " << x[j]; 
               cout << endl;    }    } 

     // Trim using the filled fragments.

     vec< vec<int> > ucores0(ucores);
     BinaryWriter::writeFile( run_dir + "/extended.long.ucores0", ucores0 );
     const int max_seg = 6;
     const double kill_ratio = 10.0;
     vec< StdMap< vec<int>, vec< triple< int, int, int > > > > segmap(max_seg+1);
     for ( int segsize = 3; segsize <= max_seg; segsize++ )
     {    vec< vec<int> > segs;
          for ( int j = 0; j < fillseqs.isize( ); j++ )
          {    const vec<int>& x = fillseqs[j];
               for ( int l = 0; l <= x.isize( ) - segsize; l++ )
               {    vec<int> y(segsize);
                    for ( int m = 1; m < segsize; m++ )
                         y[m-1] = x[l+m];
                    y[segsize-2] = x[l];
                    y[segsize-1] = x[l+segsize-1];
                    segs.push_back(y);    }    }
          ParallelSort(segs);
          for ( int i = 0; i < segs.isize( ); i++ )
          {    vec<int> key = segs[i];
               key.resize( segsize - 2 );
               int j;
               for ( j = i + 1; j < segs.isize( ); j++ )
               {    Bool key_diff = False;
                    for ( int l = 0; l < key.isize( ); l++ )
                    {    if ( segs[j][l] != key[l] )
                         {    key_diff = True;
                              break;    }    }
                    if (key_diff) break;    }
               vec< triple<int,int,int> > value;
               int max_count = 0;
               if ( j - i >= kill_ratio )
               {    for ( int m = i; m < j; m++ )
                    {    int n;
                         for ( n = m + 1; n < j; n++ )
                         {    if ( segs[n][segsize-2] != segs[m][segsize-2] ) break;
                              if ( segs[n][segsize-1] != segs[m][segsize-1] ) 
                                   break;    }
                         value.push( segs[m][segsize-2], segs[m][segsize-1], n-m );
                         max_count = Max( max_count, n-m );
                         m = n - 1;    }    }
               if ( max_count >= kill_ratio ) segmap[segsize][key] = value;
               i = j - 1;    }    }
     for ( int id = 0; id < ucores.isize( ); id++ )
          TrimByFilledFragments( ucores[id], nexts, nexts_count, segmap, True );
     BinaryWriter::writeFile( run_dir + "/extended.long.ucores", ucores );

     // Trace.

     if ( TRACE >= 0 )
     {    cout << "\nTRACE after filled fragment trimming\n" << endl;
          for ( int id = 0; id < ucores.isize( ); id++ ) 
          {    vec<int>& x = ucores[id]; 
               if ( !Member( x, TRACE ) ) continue;
               cout << "id = " << id << ", x =";
               for ( int j = 0; j < x.isize( ); j++ ) 
                    cout << " " << x[j]; 
               cout << endl;    } 
          Scram();    }

     // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

     // Load predicted gaps and use them to make a digraph.

     const int min_kmers = 800;
     cout << Date( ) << ": loading predicted gaps" << endl;
     digraphE<linklet> G;
     if (EXP2)
     {    vec< vec<int> > from(nuni), to(nuni);
          vec< vec<int> > from_edge_obj(nuni), to_edge_obj(nuni);
          vec<linklet> edges;
          String line;
          fast_ifstream in( run_dir + "/extended.unibases.k" 
               + ToString(K) + ".predicted_gaps.txt" );
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( line.Contains( "#", 0 ) ) continue;
               int u1, u2, sep, dev, nlinks;
               istringstream iline( line.c_str( ) );
               iline >> u1 >> u2 >> sep >> dev >> nlinks;
               if ( nkmers[u1] < min_kmers || nkmers[u2] < min_kmers ) continue;
               if ( u2 == to_rc[u1] ) continue; 
               /*
               if ( nlinks < MIN_LINKS_FOR_UNBRIDGEABLE )
               {    int p = BinPosition( bridge_to[u1], u2 );
                    const int max_devs = 6;
                    if ( p < 0 || bridge_to_dist[u1][p] > sep + max_devs * dev )
                    {    continue;    }    }
               */
               from[u1].push_back(u2), to[u2].push_back(u1);
               from_edge_obj[u1].push_back( edges.size( ) );
               to_edge_obj[u2].push_back( edges.size( ) );
               edges.push( sep, dev, nlinks, 0 );    }
          for ( int u = 0; u < nuni; u++ )
          {    SortSync( from[u], from_edge_obj[u] );
               SortSync( to[u], to_edge_obj[u] );    }
          G.Initialize( from, to, edges, to_edge_obj, from_edge_obj );    }

     // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

     // EXPERIMENTAL!!

     if ( EXPERIMENT || EXPERIMENT_EDIT )
     {
     cout << "\n============================================================"
          << "========================\n";
     cout << "\n" << Date( ) << ": starting experiment" << endl;
     vecbasevector longreads( data_dir + "/long_reads_orig.fastb" );

     // Define heuristics.

     const int min_flag = 100;
     const int min_flag_dist = 100;
     const int max_flag_dist = 500;

     // Find all cases.

     vec< pair<int,int> > u1_u2;
     for ( size_t id = 0; id < longreads.size( ); id++ )
     {    vec<int> xx = ucores[id];
          for ( int p1 = 0; p1 < xx.isize( ); p1++ )
          for ( int p2 = p1 + 1; p2 < xx.isize( ); p2++ )
          {    int u1 = xx[p1], u2 = xx[p2];
               if ( nkmers[u1] < min_flag || nkmers[u2] < min_flag ) continue;
               int dist = 0;
               for ( int j = p1 + 1; j < p2; j++ )
                    dist += nkmers[ xx[j] ];
               if ( dist < min_flag_dist ) continue;
               if ( dist > max_flag_dist ) continue;
               if ( Member( u1_u2, make_pair( u1, u2 ) ) ) continue;
               u1_u2.push( u1, u2 );    }    }
     PRINT( u1_u2.size( ) );

     // Execute.

     int correct_0_1 = 0, correct_1_1 = 0;
     int correct_0_2 = 0, correct_1_2 = 0, correct_2_2 = 0;
     int correct_1_3 = 0, correct_2_3 = 0, correct_3_3 = 0;
     int correct_2_4 = 0;
     int nclear = 0, nclear_wrong = 0;
     vec<String> reports( u1_u2.size( ) );
     vec< pair< vec<int>, vec<int> > > editsx;
     #pragma omp parallel for
     for ( int tt = 0; tt < u1_u2.isize( ); tt++ )
     {    if ( EXPERIMENT_TEST >= 0 && tt+1 != EXPERIMENT_TEST ) continue;
          int u1 = u1_u2[tt].first, u2 = u1_u2[tt].second;
          ostringstream fout;
          fout << "\n============================================================"
               << "========================\n\n";
          fout << Date( ) << endl;
          fout << "TEST " << tt+1 << endl;
          PRINT4_TO( fout, u1, u2, nkmers[u1], nkmers[u2] );
          fout << "\n";

          vec< triple<int,int,int> > cases;
          vec< pair<int,int> > attaches;
          vec< vec<int> > data, datax;
          for ( int id = 0; id < ucores.isize( ); id++ )
          {    const vec<int>& x = ucores[id];
               for ( int p1 = 0; p1 < x.isize( ); p1++ )
               {    if ( x[p1] != u1 ) continue;

                    const int min_attach = 10;
                    int ide = ( id < nr ? id : id - nr );
                    int attach1 = aligns[ide][ cores[id][p1] ].second.size( );
                    if ( attach1 < min_attach ) continue;

                    for ( int p2 = p1 + 1; p2 < x.isize( ); p2++ )
                    {    if ( x[p2] != u2 ) continue;

                         int attach2 = aligns[ide][ cores[id][p2] ].second.size( );
                         if ( attach2 < min_attach ) continue;

                         cases.push( id, p1, p2 );
                         attaches.push( attach1, attach2 );
                         vec<int> y;
                         y.SetToSubOf( x, p1, p2 - p1 + 1 );
                         data.push_back(y);    }    }    }
          SortSync( data, cases, attaches );
          vec< vec< triple<int,int,int> > > casesx;
          vec< vec< pair<int,int> > > attachesx;
          for ( int i = 0; i < data.isize( ); i++ )
          {    int j = data.NextDiff(i);
               datax.push_back( data[i] );
               vec< triple<int,int,int> > c;
               for ( int k = i; k < j; k++ )
                    c.push_back( cases[k] );
               casesx.push_back(c);
               vec< pair<int,int> > a;
               for ( int k = i; k < j; k++ )
                    a.push_back( attaches[k] );
               attachesx.push_back(a);
               i = j - 1;    }
          int nx = datax.size( );
          PRINT2_TO( fout, data.size( ), nx );
          fout << "\n";
          vec<basevector> datab, datab0;
          for ( int i = 0; i < nx; i++ )
          {    const vec<int>& x = datax[i];
               basevector b = unibases[ x[0] ];
               for ( int j = 1; j < x.isize( ); j++ )
               {    b.resize( b.isize( ) - (K-1) );
                    b = Cat( b, unibases[ x[j] ] );    }
               basevector b0;
               if ( x.isize( ) >= 2 )
               {    b0 = unibases[ x[1] ];
                    for ( int j = 2; j < x.isize( ) - 1; j++ )
                    {    b0.resize( b0.isize( ) - (K-1) );
                         b0 = Cat( b0, unibases[ x[j] ] );    }    }
               datab.push_back(b), datab0.push_back(b0);    }
          vec<Bool> valid( nx, False );
          for ( int j = 0; j < nx; j++ )
          {    fout << "[" << j+1 << "," << casesx[j].size( ) << "x]";
               const vec<int>& x = datax[j];
               for ( int l = 0; l < x.isize( ); l++ )
                    fout << " " << x[l];
               if (VALIDATE)
               {    vec<placementy> places = FindGenomicPlacementsY( 
                         0, datab[j], LG, genome2, Glocs );
                    if ( places.nonempty( ) ) 
                    {    valid[j] = True;
                         fout << " [TRUE]";    }    }
               fout << "\n";    }

          // Check for very low coverage.

          const int min_count = 2;
          int best = 0;
          for ( int j = 0; j < nx; j++ )
               best = Max( best, casesx[j].isize( ) );
          if ( best < min_count )
          {    if (VALIDATE)
               {    if ( Sum(valid) == 0 && nx == 1 )
                    {    
                         #pragma omp critical
                         {    correct_0_1++;    }
                         continue;    }
                    if ( Sum(valid) == 0 && nx == 2 )
                    {    
                         #pragma omp critical
                         {    correct_0_2++;    }
                         continue;    }    }
               fout << "Coverage too low, treating as junk." << endl;
               reports[tt] = fout.str( );
               continue;    }

          if ( datax.solo( ) )
          {    fout << "There is only one answer!\n";
               if ( VALIDATE && valid[0] ) 
               {
                    #pragma omp critical
                    {    correct_1_1++;    }    }
               else reports[tt] = fout.str( );
               continue;    }
          ostringstream dout;
          vec< vec< vec< pair<ho_interval,ho_interval> > > > signal(nx);
          for ( int j = 0; j < nx; j++ )
               signal[j].resize(nx);
          for ( int j1 = 0; j1 < nx; j1++ )
          for ( int j2 = j1+1; j2 < nx; j2++ )
          {    dout << "\n[" << j1+1 << "] vs [" << j2+1 << "]\n";
               if ( datab0[j1].size( ) == 0 || datab0[j2].size( ) == 0 )
                    dout << "\n(indel)\n";
               else
               {    align a;
                    SmithWatFreeSym( datab0[j1], datab0[j2], a, True, True, 1, 1 );
                    PrintVisualAlignment( True, dout, datab0[j1], datab0[j2], a ); 

                    // Find 'signal' windows defined by an isolated mismatch.
                    // (generalized)
     
                    const int flank = 5;
                    const int max_hop = 5;
                    vec<ho_interval> perfs1, perfs2;
                    a.PerfectIntervals1( datab0[j1], datab0[j2], perfs1 );
                    a.PerfectIntervals2( datab0[j1], datab0[j2], perfs2 );
                    for ( int l = 0; l < perfs1.isize( ) - 1; l++ )
                    {    if ( perfs1[l].Length( ) < flank ) continue;
                         for ( int m = l + 1; m < perfs1.isize( ); m++ )
                         {    if ( perfs1[m].Length( ) < flank ) continue;
                              if ( perfs1[m].Start( ) - perfs1[l].Stop( ) > max_hop )
                                   continue;
                              if ( perfs2[m].Start( ) - perfs2[l].Stop( ) > max_hop )
                                   continue;
                              ho_interval h1( perfs1[l].Stop( ) - flank, 
                                   perfs1[m].Start( ) + flank );
                              ho_interval h2( perfs2[l].Stop( ) - flank, 
                                   perfs2[m].Start( ) + flank );
                              dout << "signal: " << h1 << " vs " << h2 << "\n";    
                              signal[j1][j2].push(h1,h2), signal[j2][j1].push(h2,h1);
                              break;    }    }    }    }
          vec< vec< pair<int,int> > > cindex0(nuni);
          for ( int id = 0; id < ucores.isize( ); id++ )
          for ( int p = 0; p < ucores[id].isize( ); p++ )
               cindex0[ ucores[id][p] ].push( id, p );
          vec<int> score( nx, 0 ), wins( nx, 0 ), wins_weighted( nx, 0 );
          vec< vec<int> > nots_all;
          for ( int j1 = 0; j1 < nx; j1++ )
          {    const int infinity = 1000000000;
               for ( int l1 = 0; l1 < casesx[j1].isize( ); l1++ )
               {    int id = casesx[j1][l1].first;
                    int p1 = casesx[j1][l1].second, p2 = casesx[j1][l1].third;
                    align a;
                    Bool rc = ( id >= nr );
                    int ide = ( rc ? id - nr : id ), rpos1, rpos2, upos1, upos2;
                    const vec< pair<int,int> >& a1 
                         = aligns[ide][ cores[id][p1] ].second;
                    const vec< pair<int,int> >& a2 
                         = aligns[ide][ cores[id][p2] ].second;
                    if ( !rc )
                    {    int ptr1 = a1.isize( ) - 1, ptr2 = 0;
                         while(1)
                         {    int rpos1 = a1[ptr1].first, rpos2 = a2[ptr2].first;
                              if ( rpos2 >= rpos1 ) break;
                              if ( ptr1 == 0 && ptr2 == a2.isize( ) - 1 ) break;
                              if ( ptr1 > 0 ) ptr1--;
                              if ( ptr2 < a2.isize( ) - 1 ) ptr2++;    }
                         rpos1 = a1[ptr1].first, upos1 = a1[ptr1].second;
                         rpos2 = a2[ptr2].first, upos2 = a2[ptr2].second;    }
                    else
                    {    int ptr1 = a2.isize( ) - 1, ptr2 = 0;
                         while(1)
                         {    int rpos1 = a2[ptr1].first, rpos2 = a1[ptr2].first;
                              if ( rpos2 >= rpos1 ) break;
                              if ( ptr1 == 0 && ptr2 == a1.isize( ) - 1 ) break;
                              if ( ptr1 > 0 ) ptr1--;
                              if ( ptr2 < a1.isize( ) - 1 ) ptr2++;    }
                         rpos1 = a2[ptr1].first, upos1 = a2[ptr1].second;
                         rpos2 = a1[ptr2].first, upos2 = a1[ptr2].second;    }

                    vec<int> errsx( nx, 0 );
                    int rsize = rpos2 - rpos1;
                    vec< pair<int,int> > sig;
                    for ( int j = 0; j < nx; j++ )
                    {    if ( !rc )
                         {    int zsize = datab[j].isize( ) - upos1
                                   - ( unibases[u2].isize( ) - upos2 );
                              if ( rsize < 0 ) errsx[j] = zsize - rsize;
                              else if ( zsize < 0 ) errsx[j] = rsize - zsize;
                              else if ( rsize == 0 || zsize == 0 )
                              {    errsx[j] = rsize + zsize;    }
                              else
                              {    basevector z, r = longreads[ide];
                                   r.SetToSubOf( r, rpos1, rsize );
                                   z.SetToSubOf( datab[j], upos1, zsize );
                                   errsx[j] = SmithWatFreeSym(
                                        r, z, a, True, True, 1, 1 );    
                                   vec<ho_interval> perfs2;
                                   a.PerfectIntervals2( r, z, perfs2 );
                                   for ( int l = 0; l < perfs2.isize( ); l++ )
                                   {    perfs2[l] = perfs2[l] + upos1
                                             - nkmers[ datax[j][0] ];    }
                                   for ( int l = 0; l < nx; l++ )
                                   {    for ( int m = 0; m < signal[j][l].isize( ); 
                                             m++ )
                                        {    ho_interval h = signal[j][l][m].first;
                                             for ( int n = 0; n < perfs2.isize( ); 
                                                  n++ )
                                             {    if ( Subset( h, perfs2[n] ) )
                                                       sig.push( j, l );
                                        }    }    }    }    }

                    // ---->

                    else
                    {    int zsize = datab[j].isize( ) - upos1
                              - ( unibases[u1].isize( ) - upos2 );
                         if ( rsize < 0 ) errsx[j] = zsize - rsize;
                         else if ( zsize < 0 ) errsx[j] = rsize - zsize;
                         else if ( rsize == 0 || zsize == 0 )
                         {    errsx[j] = rsize + zsize;    }
                         else
                         {    basevector z, r = longreads[ide], d = datab[j];
                              r.SetToSubOf( r, rpos1, rsize );
                              d.ReverseComplement( );
                              z.SetToSubOf( d, upos1, zsize );
                              errsx[j] = SmithWatFreeSym(
                                   r, z, a, True, True, 1, 1 );    

                              vec<ho_interval> perfs2;
                              a.PerfectIntervals2( r, z, perfs2 );
                              for ( int l = 0; l < perfs2.isize( ); l++ )
                              {    int start = perfs2[l].Start( );
                                   int stop = perfs2[l].Stop( );
                                   int xstart = datab[j].isize( ) - upos1 - stop
                                        - nkmers[ datax[j][0] ];
                                   int xstop = datab[j].isize( ) - upos1 - start
                                        - nkmers[ datax[j][0] ];
                                   perfs2[l] = ho_interval( 
                                        xstart, xstop );    }

                              for ( int l = 0; l < nx; l++ )
                              {    for ( int m = 0; m < signal[j][l].isize( ); 
                                        m++ )
                                   {    ho_interval h = signal[j][l][m].first;
                                        for ( int n = 0; n < perfs2.isize( ); 
                                             n++ )
                                        {    if ( Subset( h, perfs2[n] ) )
                                                  sig.push( j, l );
                                   }    }    }    }    }    }
               dout << "\ntesting [" << j1+1 << "." << l1+1 << "=" << id 
                    << "] -- " << "see";
               for ( int j = 0; j < nx; j++ )
               {    dout << " ";
                    if ( errsx[j] == infinity ) dout << "inf";
                    else dout << errsx[j];
                    score[j] += errsx[j];    }
               int m = Min(errsx);
               vec<int> winners;
               for ( int j = 0; j < nx; j++ )
                    if ( errsx[j] == m ) winners.push_back(j);
               if ( winners.solo( ) ) 
               {    wins[ winners[0] ]++;
                    pair<int,int> a = attachesx[j1][l1];
                    wins_weighted[ winners[0] ] 
                         += Min( a.first, a.second );    }
               dout << "\n";    
               vec<int> nots;
               for ( int j = 0; j < sig.isize( ); j++ )
                    nots.push_back( sig[j].second );
               UniqueSort(nots);
               if ( nots.nonempty( ) )
               {    nots_all.push_back(nots);
                    dout << "not ";
                    for ( int j = 0; j < nots.isize( ); j++ )
                    {    if ( j > 0 ) dout << " or ";
                         dout << nots[j] + 1;    }
                    dout << "\n";    }    }    }

     // Using the weights, declare certain losers, and adjust nots_all 
     // accordingly.

     const double min_ratio = 20.0;
     int mw = Max(wins_weighted);
     vec<Bool> loser( nx, False );
     for ( int j = 0; j < nx; j++ )
          if ( wins_weighted[j] * min_ratio < mw ) loser[j] = True;
     Sort(nots_all);
     for ( int j = 0; j < nots_all.isize( ); j++ )
     {    for ( int l = 0; l < nx; l++ )
          {    if ( loser[l] && !Member( nots_all[j], l ) )
                    nots_all[j].push_back(l);    }    }

     // Summary analysis.

     int win_count = 0;
     for ( int j = 0; j < nx; j++ )
          if ( wins[j] > 0 ) win_count++;
     vec<int> strong_votes( nx, 0 );
     for ( int j = 0; j < nots_all.isize( ); j++ )
     {    if ( nots_all[j].isize( ) == nx - 1 )
          {    for ( int k = 0; k < nx; k++ )
               {    if ( !Member( nots_all[j], k ) )
                         strong_votes[k]++;    }    }    }
     Bool power = True;
     for ( int j1 = 0; j1 < nx; j1++ )
     for ( int j2 = j1+1; j2 < nx; j2++ )
          if ( signal[j1][j2].empty( ) ) power = False;
     const int min_sig = 4;
     Bool sig_solved = False;
     int strong_voters = 0;
     if (power)
     {    sig_solved = True;
          for ( int j = 0; j < nx; j++ )
          {    if ( strong_votes[j] >= min_sig ) strong_voters++;
               if ( strong_votes[j] > 0 && strong_votes[j] < min_sig )
                    sig_solved = False;    }    }

     // Identify a class of clear-cut cases.  Note that there are some false
     // positives.  All known false positives are cases where we would convert
     // one wrong sequence into another wrong sequence.

     Bool clear_cut = False;
     if ( strong_voters == 1 && sig_solved ) clear_cut = True;
     if (clear_cut)
     {    
          #pragma omp critical
          {    if ( strong_votes[0] > 0 ) editsx.push( datax[1], datax[0] );
               else editsx.push( datax[0], datax[1] );
               nclear++;    }
          if (VALIDATE)
          {    Bool OK = True;
               for ( int j = 0; j < nx; j++ )
               {    if ( strong_votes[j] > 0 && !valid[j] ) OK = False;
                    if ( strong_votes[j] == 0 && valid[j] ) OK = False;    }
               if ( !OK ) 
               {    fout << "\nOops, we're making a mistake.\n";
                    #pragma omp critical
                    {    nclear_wrong++;    }    }    }    }

     // Other cases.

     if ( VALIDATE && sig_solved )
     {    Bool OK = True;
          for ( int j = 0; j < nx; j++ )
          {    if ( strong_votes[j] > 0 && !valid[j] ) OK = False;
               if ( strong_votes[j] == 0 && valid[j] ) OK = False;    }
          if (OK)
          {    
               if ( strong_voters == 1 && nx == 2 )
               {    
                    #pragma omp critical
                    {    correct_1_2++;    }
                    continue;    }
               if ( strong_voters == 2 && nx == 2 )
               {    
                    #pragma omp critical
                    {    correct_2_2++;    }
                    continue;    }
               if ( strong_voters == 1 && nx == 3 )
               {    
                    #pragma omp critical
                    {    correct_1_3++;    }
                    continue;    }
               if ( strong_voters == 2 && nx == 3 )
               {    
                    #pragma omp critical
                    {    correct_2_3++;    }
                    continue;    }
               if ( strong_voters == 3 && nx == 3 )
               {    
                    #pragma omp critical
                    {    correct_3_3++;    }
                    continue;    }    }    }

     Bool weighted_win = True;
     int mwi = -1;
     for ( int j = 0; j < nx; j++ )
          if ( wins_weighted[j] == mw ) mwi = j;
     for ( int j = 0; j < nx; j++ )
     {    if ( j == mwi ) continue;
          if ( wins_weighted[j] * min_ratio >= mw ) weighted_win = False;    }
     if ( VALIDATE && weighted_win )
     {    Bool OK = True;
          for ( int j = 0; j < nx; j++ )
          {    if ( j == mwi && !valid[j] ) OK = False;
               if ( j != mwi && valid[j] ) OK = False;    }
          if ( OK && nx == 2 )
          {    
               #pragma omp critical
               {    correct_1_2++;    }
               continue;    }
          if ( OK && nx == 3 )
          {    
               #pragma omp critical
               {    correct_1_3++;    }
               continue;    }    }

     Bool solved = False;
     if ( win_count == 1 || sig_solved ) solved = True;

     if ( VALIDATE && win_count == 1 && nx == 2 )
     {    Bool OK = True;
          for ( int j = 0; j < nx; j++ )
          {    if ( wins[j] > 0 && !valid[j] ) OK = False;
               if ( wins[j] == 0 && valid[j] ) OK = False;    }
          if (OK) 
          {    
               #pragma omp critical
               {    correct_1_2++;    }
               continue;    }    }

     ostringstream xout;
     xout << fout.str( );
     if ( !solved || EXPERIMENT_VERBOSITY >= 2 ) xout << dout.str( );

     if ( win_count == 1 ) xout << "Unanimous!\n";

     xout << "\nSCORE SUMMARY: ";
     score.Print(xout);
     xout << "\nWINNERS:       ";
     wins.Print(xout);
     xout << "\nWEIGHTED:      ";
     wins_weighted.Print(xout);
     xout << "\nSTRONG VOTES:  ";
     strong_votes.Print(xout);
     xout << "\n";
     if (sig_solved) xout << "\nsignal answers question\n";
     reports[tt] = xout.str( );

     }

     // Print reports.

     if ( EXPERIMENT_VERBOSITY >= 1 )
     {    for ( int j = 0; j < reports.isize( ); j++ )
               cout << reports[j];    }

     // Print summary.

     cout << "\n============================================================"
          << "========================\n\n";
     if (VALIDATE)
     {    cout << "GRAND SUMMARY\n";
          cout << correct_0_1 << " cases where we correctly picked 0 of 1" << endl;
          cout << correct_1_1 << " cases where we correctly picked 1 of 1" << endl;
          cout << correct_0_2 << " cases where we correctly picked 0 of 2" << endl;
          cout << correct_1_2 << " cases where we correctly picked 1 of 2" << endl;
          cout << correct_2_2 << " cases where we correctly picked 2 of 2" << endl;
          cout << correct_1_3 << " cases where we correctly picked 1 of 3" << endl;
          cout << correct_2_3 << " cases where we correctly picked 2 of 3" << endl;
          cout << correct_3_3 << " cases where we correctly picked 3 of 3" << endl;
          cout << u1_u2.isize( ) - correct_0_1 - correct_1_1 - correct_0_2 
               - correct_1_2 - correct_2_2 - correct_1_3 - correct_2_3 - correct_3_3
               << " other cases\n\n";    }
     if (VALIDATE) 
     {    PRINT2( nclear, nclear_wrong );    }
     else PRINT(nclear);
     cout << "\n" << Date( ) << ": experiment complete" << endl;
     if ( !EXPERIMENT_EDIT ) Scram();

     // Make edits. 

     cout << Date( ) << ": making experimental edits" << endl;
     int edit_count = 0;
     for ( int id = 0; id < ucores.isize( ); id++ )
     {    vec<int>& x = ucores[id];
          for ( int j = 0; j < editsx.isize( ); j++ )
          {    for ( int p = 0; p < x.isize( ); p++ )
               {    if ( x.Contains( editsx[j].first, p ) )
                    {    edit_count++;
                         vec<int> x1, x2, x3;
                         x1.SetToSubOf( x, 0, p );
                         x2 = editsx[j].second;
                         x3.SetToSubOf( x, p + editsx[j].first.isize( ),
                              x.isize( ) - p - editsx[j].first.isize( ) );
                         x = x1;
                         x.append(x2);
                         x.append(x3);    
                         break;    }    }    }    }
     PRINT(edit_count);

     }

     // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

     // Compute cindex.

     cout << Date( ) << ": compute cindex" << endl;
     vec< vec< pair<int,int> > > cindex(nuni);
     for ( int id = 0; id < ucores.isize( ); id++ )
     for ( int p = 0; p < ucores[id].isize( ); p++ )
          cindex[ ucores[id][p] ].push( id, p );

     // Set up useq class.

     useq dummy;
     dummy.SetUnibases( K, unibases, to_rc );

     // Build neighborhoods.

     const int max_seed = 500;
     vec< vec<int> > accepted;
     vec<String> reports(nuni);
     Bool verbose = False;
     const int min_radius = 300;
     int nwarnings = 0;
     vec< vec< vec<int> > > lexts(nuni), rexts(nuni);
     #pragma omp parallel for
     for ( int test_id = 0; test_id < nuni; test_id++ )
     {    if ( TEST_ID >= 0 && test_id != TEST_ID ) continue;
          if ( nkmers[test_id] < MIN_SEED ) continue;
          ostringstream hout;
          for ( int pass = 1; pass <= 2; pass++ )
          // for ( int pass = 1; pass <= 3; pass++ )
          {    // if ( pass == 2 && nkmers[test_id] <= max_seed ) continue;
     
               if ( pass == 3 )
               {    hout << "\n";
                    vec< pair< vec<Bool>, vec<Bool> > > vals;
                    for ( int v = 0; v < cindex[test_id].isize( ); v++ )
                    {    int id1 = cindex[test_id][v].first;
                         const vec<int>& x = ucores[id1];
                         int p = cindex[test_id][v].second;
                         vec<Bool> lm( lexts[test_id].size( ), False );
                         vec<Bool> rm( rexts[test_id].size( ), False );
                         for ( int j = 0; j < lexts[test_id].isize( ); j++ )
                         {    if ( Match( lexts[test_id][j], 
                                   lexts[test_id][j].isize( ) - 1, x, p ) )
                              {    lm[j] = True;    }    }
                         for ( int j = 0; j < rexts[test_id].isize( ); j++ )
                         {    if ( Match( rexts[test_id][j], 0, x, p ) )
                              {    rm[j] = True;    }    }
                         if ( Sum(lm) == 0 || Sum(rm) == 0 ) continue;
                         vals.push( lm, rm );    }
                    Sort(vals);
                    for ( int i = 0; i < vals.isize( ); i++ )
                    {    int j = vals.NextDiff(i);
                         const vec<Bool> &lm = vals[i].first, &rm = vals[i].second;
                         hout << "count = " << j-i << "; left exts =";
                         for ( int j = 0; j < lm.isize( ); j++ )
                              if ( lm[j] ) hout << " " << j+1;
                         hout << "; right exts =";
                         for ( int j = 0; j < rm.isize( ); j++ )
                              if ( rm[j] ) hout << " " << j+1;
                         hout << "\n";
                         i = j - 1;    }
                    break;    }

               for ( int radius = MAX_RADIUS; radius >= min_radius; radius -= 100 )
               {    int left_flank = radius, right_flank = radius;
                    // if ( nkmers[test_id] > max_seed )
                    {    if ( pass == 1 ) left_flank = 0;
                         else right_flank = 0;    }
                    vec< pair< vec<int>, int > > hypos;
                    vec<int> scores, valids;
                    for ( int m = 0; m < cindex[test_id].isize( ); m++ )
                    {    int id1 = cindex[test_id][m].first;
                         int p = cindex[test_id][m].second;
                         const vec<int>& x = ucores[id1];
                         int left_kmers = 0, right_kmers = 0, j1, j2;
                         if ( left_flank == 0 ) j1 = p;
                         else
                         {    for ( j1 = p - 1; j1 >= 0; j1-- )
                              {    left_kmers += nkmers[ x[j1] ];
                                   if ( left_kmers >= left_flank ) break;    }    }
                         if ( right_flank == 0 ) j2 = p;
                         else
                         {    for ( j2 = p + 1; j2 < x.isize( ); j2++ )
                              {    right_kmers += nkmers[ x[j2] ];
                                   if ( right_kmers >= right_flank ) break;    }    }
                         if ( left_kmers < left_flank || right_kmers < right_flank ) 
                              continue;
                         vec<int> h;
                         h.SetToSubOf( x, j1, j2 - j1 + 1 );
                         hypos.push( h, p-j1 );    }
                    Sort(hypos);

                    vec<Bool> hypos_to_delete( hypos.size( ), False );
                    for ( int i = 0; i < hypos.isize( ); i++ )
                    {    int j;
                         for ( j = i + 1; j < hypos.isize( ); j++ )
                              if ( hypos[j].first != hypos[i].first ) break;
                         const vec<int>& x = hypos[i].first;
                         Bool bad = False;
                         for ( int l = 1; l < illegal.isize( ); l++ )
                         {    if (bad) break;
                              for ( int j = 0; j < illegal[l].isize( ); j++ )
                              {    if ( x.Contains( illegal[l][j] ) )
                                   {    bad = True;
                                        break;    }    }    }
                         if (bad)
                         {    for ( int k = i; k < j; k++ )
                                   hypos_to_delete[k] = True;    }
                         i = j - 1;    }
                    EraseIf( hypos, hypos_to_delete );

                    vec< vec<int> > hypos2;
                    vec<int> hypo_pos, hypocounts;
                    for ( int i = 0; i < hypos.isize( ); i++ )
                    {    int j = hypos.NextDiff(i);
                         hypocounts.push_back(j-i);
                         hypos2.push_back( hypos[i].first );
                         hypo_pos.push_back( hypos[i].second );
                         i = j - 1;    }
                    for ( int j = 0; j < hypos2.isize( ); j++ )
                    {    const vec<int>& x1 = hypos2[j];
                         basevector b = unibases[ x1[0] ];
                         for ( int z = 1; z < x1.isize( ); z++ )
                         {    b.resize( b.isize( ) - (K-1) );
                              b = Cat( b, unibases[ x1[z] ] );    }
                         int valid = 0;
                         if (VALIDATE)
                         {    vec<placementy> places = FindGenomicPlacementsY( 
                                   0, b, LG, genome2, Glocs );
                              for ( int l = 0; l < places.isize( ); l++ )
                              {    if ( places[l].pos 
                                        < genome[ places[l].g ].isize( ) ) 
                                   {    valid++;    }    }    }
                         int n1 = 0;
                         for ( int z = 0; z < x1.isize( ); z++ )
                              n1 += nkmers[ x1[z] ];
                         if (verbose)
                         {    hout << "\n" << String( 100, '=' ) << "\n\n";
                              PRINT_TO( hout, n1 );
                              hout << "x1 =";
                              for ( int j = 0; j < x1.isize( ); j++ )
                                   hout << " " << x1[j];
                              hout << "\n";    }
                         vec< pair<ualign, int> > aligns_ids;
                         int id1 = 0; // *******************************************
	                 GetAligns( id1, hypo_pos[j], x1, nkmers, ucores, 
                              cindex, aligns_ids );
                         int test_cov = 0;
                         for ( int z = 0; z < aligns_ids.isize( ); z++ )
                         {    const ualign& a = aligns_ids[z].first;
                              int id2 = aligns_ids[z].second;
                              const vec<int>& x2 = ucores[id2];
                              useq u1(x1), u2(x2);
                              if (verbose)
                              {    hout << "\n";
                                   PRINT_TO( hout, id2 );
                                   a.Print( hout, u1, u2 );    }
                              vec<Bool> t( x1.size( ), False );
                              for ( int z = 0; z < a.Ties( ).isize( ); z++ )
                                   t[ a.Tie(z).first ] = True;
                              if ( t[ hypo_pos[j] ]
                                   && ( left_flank == 0 || t[ hypo_pos[j] - 1 ] )
                                   && ( right_flank == 0 || t[ hypo_pos[j] + 1 ] ) )
                              {    test_cov++;    }    }
                         scores.push_back(test_cov * x1.isize( ));
                         valids.push_back(valid);
                         if (verbose)
                         {    PRINT3_TO( hout,
                                   valid, test_cov, hypocounts[j] );    }    }
                    ReverseSortSync( scores, hypos2, hypo_pos, hypocounts, valids );
                    const int min_score = 5;
                    if ( hypocounts.empty( ) || Max(hypocounts) < min_score ) 
                    {    if ( radius == min_radius )
                         {    hout << "\ntesting unibase " << test_id << ", " 
                                   << nkmers[test_id] << " kmers, " 
                                   << uperfect[test_id] << "x" << ", radius = " 
                                   << radius;
                              if ( left_flank == 0 ) hout << ", looking right";
                              if ( right_flank == 0 ) hout << ", looking left";
                              hout << "\n";
                              PRINT_TO( hout, hypocounts.size( ) );
                              if ( hypocounts.nonempty( ) )
                                   PRINT2_TO( hout, Max(hypocounts), min_score );
                              if (VALIDATE)
                              {    int gcount = TrueLocs( test_id, left_flank, 
                                        right_flank, K, unibases, nkmers, nexts, 
                                        to_rc, genome, genome2, LG, Glocs );
                                   if ( gcount > 0 )
                                   {    hout << "WARNING: Failed to find true site." 
                                             << endl; 
                                        #pragma omp critical
                                        {    nwarnings++;    }
                                        PRINT_TO( hout, gcount );    }    }    }
                         continue;    }
     
                    int gcount = 0;
                    if (VALIDATE)
                    {    gcount = TrueLocs( test_id, left_flank, right_flank, K, 
                              unibases, nkmers, nexts, to_rc, genome, genome2, 
                              LG, Glocs );   }
          
                    if (verbose)
                    {    hout << "\n#####################################"
                              << "##########################\n\n";    }
                    else 
                    {    hout << "\ntesting unibase " << test_id << ", " 
                              << nkmers[test_id] << " kmers, " << uperfect[test_id] 
                              << "x" << ", radius = " << radius;
                         if ( left_flank == 0 ) hout << ", looking right";
                         if ( right_flank == 0 ) hout << ", looking left";
                         hout << endl;    }
                    int hcount = 0, finds = 0;
                    int min_count_valid = 1000000;

                    int best_cscore = 0, best_cscore0 = 0;
                    for ( int j = 0; j < hypos2.isize( ); j++ )
                    {    if ( scores[j] * MC_REL < scores[0] ) break;
                         best_cscore = Max( best_cscore, hypocounts[j] 
                              * hypos2[j].isize( ) * hypos2[j].isize( ) );
                         best_cscore0 = Max( best_cscore0, hypocounts[j] 
                              * hypos2[j].isize( ) );    }

                    for ( int j = 0; j < hypos2.isize( ); j++ )
                    {    if ( scores[j] * MC_REL < scores[0] ) break;

                         int cscore =  hypocounts[j] 
                              * hypos2[j].isize( ) * hypos2[j].isize( );
                         int cscore0 =  hypocounts[j] * hypos2[j].isize( );
                         if ( ( cscore * MC_BOUND < best_cscore )
                              && ( cscore0 * MC_BOUND < best_cscore0 ) )
                         {    continue;    }

                         hout << "[" << ++hcount << "] ";
                         const vec<int>& x1 = hypos2[j];
                         if ( pass == 1 ) rexts[test_id].push_back(x1);
                         if ( pass == 2 ) lexts[test_id].push_back(x1);
                         #pragma omp critical
                         {    accepted.push_back(x1);    }
                         hout << "score = " << scores[j] << ", count = " 
                              << hypocounts[j] << ", valid = " << valids[j] 
                              << ", x1 =";
                         if ( valids[j] > 0 ) 
                              min_count_valid = Min(min_count_valid, hypocounts[j]);
                         finds += valids[j];
                         for ( int z = 0; z < x1.isize( ); z++ )
                              hout << " " << x1[z];
                         hout << "\n";    }

                    if ( hypocounts.nonempty( ) )
                    {    PRINT3_TO( hout, Max(hypocounts), min_count_valid,
                              double( Max(hypocounts) ) 
                              / double(min_count_valid) );    }

                    // Test number of finds.
          
                    if (VALIDATE)
                    {    if ( finds > gcount )
                         {    hout << "ERROR: finds > gcount" << endl;
                              hout << "This is probably because we're at the end "
                                   << "of a chromosome, not a real problem." << endl;
                              PRINT2_TO( hout, finds, gcount );    }
                         if ( finds < gcount )
                         {    hout << "WARNING: Failed to find true site." << endl; 
                              #pragma omp critical
                              {    nwarnings++;    }
                              PRINT2_TO( hout, finds, gcount );    }    }
                    break;    }    }
          reports[test_id] = hout.str( );    }

     // Print reports.

     if (VERBOSE)
     {    for ( int test_id = 0; test_id < nuni; test_id++ )
               cout << reports[test_id];    }

     // Hard filter.

     cout << Date( ) << ": start hard filter" << endl;
     while(1)
     {    int dels = 0;
          for ( int u1 = 0; u1 < nuni; u1++ )
          {    vec<Bool> ldel( lexts[u1].size( ), False );
               for ( int j1 = 0; j1 < lexts[u1].isize( ); j1++ )
               {    const vec<int>& x1 = lexts[u1][j1];
                    for ( int p1 = 1; p1 < x1.isize( ); p1++ )
                    {    int u2 = x1[p1];
                         if ( nkmers[u2] < MIN_SEED ) continue;
                         Bool match = False;
                         for ( int j2 = 0; j2 < rexts[u2].isize( ); j2++ )
                         {    if ( Match( x1, p1, rexts[u2][j2], 0 ) )
                              {    match = True;
                                   break;    }    }
                         if ( !match ) ldel[j1] = True;
                         match = False;
                         for ( int j2 = 0; j2 < lexts[u2].isize( ); j2++ )
                         {    if ( Match( x1, p1, lexts[u2][j2], 
                                   lexts[u2][j2].isize( ) - 1 ) )
                              {    match = True;
                                   break;    }    }
                         if ( !match ) ldel[j1] = True;    }    }
               dels += Sum(ldel);
               EraseIf( lexts[u1], ldel );
               vec<Bool> rdel( rexts[u1].size( ), False );
               for ( int j1 = 0; j1 < rexts[u1].isize( ); j1++ )
               {    const vec<int>& x1 = rexts[u1][j1];
                    for ( int p1 = 0; p1 < x1.isize( ) - 1; p1++ )
                    {    int u2 = x1[p1];
                         if ( nkmers[u2] < MIN_SEED ) continue;
                         Bool match = False;
                         for ( int j2 = 0; j2 < rexts[u2].isize( ); j2++ )
                         {    if ( Match( x1, p1, rexts[u2][j2], 0 ) )
                              {    match = True;
                                   break;    }    }
                         if ( !match ) rdel[j1] = True;
                         match = False;
                         for ( int j2 = 0; j2 < lexts[u2].isize( ); j2++ )
                         {    if ( Match( x1, p1, lexts[u2][j2], 
                                   lexts[u2][j2].isize( ) - 1 ) )
                              {    match = True;
                                   break;    }    }
                         if ( !match ) rdel[j1] = True;    }    }
               dels += Sum(rdel);
               EraseIf( rexts[u1], rdel );    }
          PRINT(dels);
          if ( dels == 0 ) break;    }

     // Now we can extend each right neighborhood by the common part of the left
     // neighborhoods, and conversely.

     accepted.clear( );
     /*
     for ( int u = 0; u < nuni; u++ )
     {    int lshare = 1, rshare = 1;
          for ( int j = 0; j < lexts[u].isize( ); j++ )
          {    if ( lexts[u][j].isize( ) < lshare + 1 ) break;
               if ( lexts[u][j][ lexts[u][j].isize( ) - lshare - 1 ]
                    != lexts[u][0][ lexts[u][0].isize( ) - lshare - 1 ] )
               {    break;    }
               lshare++;    }
          for ( int j = 0; j < rexts[u].isize( ); j++ )
          {    if ( rexts[u][j].isize( ) < rshare + 1 ) break;
               if ( rexts[u][j][rshare] != rexts[u][0][rshare] ) break;
               rshare++;    }
          for ( int j = 0; j < lexts[u].isize( ); j++ )
          {    vec<int> x = lexts[u][j];
               for ( int m = 1; m < rshare; m++ )
                    x.push_back( rexts[u][0][m] );
               if (VERBOSE)
               {    cout << "saving";
                    for ( int j = 0; j < x.isize( ); j++ )
                         cout << " " << x[j];
                    cout << "\n";    }
               accepted.push_back(x);    }
          for ( int j = 0; j < rexts[u].isize( ); j++ )
          {    vec<int> x = rexts[u][j];
               for ( int m = 1; m < lshare; m++ )
                    x.push_front( lexts[u][0][ lexts[u][0].isize( ) - m - 1 ] );
               if (VERBOSE)
               {    cout << "saving";
                    for ( int j = 0; j < x.isize( ); j++ )
                         cout << " " << x[j];
                    cout << "\n";    }
               accepted.push_back(x);    }    }
     */
     for ( int u = 0; u < nuni; u++ )
     {    accepted.append( lexts[u] );
          accepted.append( rexts[u] );    }

     // NEW EXPERIMENT!!

     vec< vec<int> > exp2_verts;
     if (EXP2)
     {

     // Build a graph from the accepteds.

     vec< vec<int> > acc(accepted);
     for ( int i = 0; i < accepted.isize( ); i++ )
     {    vec<int> x = accepted[i];
          x.ReverseMe( );
          for ( int j = 0; j < x.isize( ); j++ )
               x[j] = to_rc[ x[j] ];
          acc.push_back(x);    }
     UniqueSort(acc);
     vec< vec< pair<int,int> > > aindex(nuni);
     for ( int id = 0; id < acc.isize( ); id++ )
     for ( int p = 0; p < acc[id].isize( ); p++ )
          aindex[ acc[id][p] ].push( id, p );
     vec< vec<int> > from( acc.size( ) ), to( acc.size( ) );
     for ( int id1 = 0; id1 < acc.isize( ); id1++ )
     {    const vec<int>& x = acc[id1];
          for ( int j = 0; j < aindex[ x.back( ) ].isize( ); j++ )
          {    int id2 = aindex[ x.back( ) ][j].first;
               int p = aindex[ x.back( ) ][j].second;
               const vec<int>& y = acc[id2];
               if ( p == y.isize( ) - 1 ) continue;
               Bool mismatch = False;
               int overlap = 0;
               for ( int l = Max( 0, p - ( x.isize( ) - 1 ) ); l <= p; l++ )
               {    if ( x[ l - p + ( x.isize( ) - 1 ) ] != y[l] )
                    {    mismatch = True;
                         break;    }
                    overlap += nkmers[ y[l] ];    }
               if (mismatch) continue;
               if ( overlap < KBIG-K ) continue;
               from[id1].push_back(id2), to[id2].push_back(id1);    }    }
     for ( int i = 0; i < acc.isize( ); i++ )
     {    Sort(from[i]), Sort(to[i]);    }
     digraph A( from, to );

     // Now start EXP2.

     const int stud_size = 1000;
     ForceAssertGe( stud_size, KBIG );
     vec< vec< pair<int,int> > > cindex1(nuni);
     for ( int id = 0; id < ucores0.isize( ); id++ )
     for ( int p = 0; p < ucores0[id].isize( ); p++ )
          cindex1[ ucores0[id][p] ].push( id, p );

     // Go through the unipath studs and walk right to the next stud.
     // In the process we build "units".  We define a unit to be a sequence
     // x1,...xm,y1,...,yn of unipaths such that 
     //     nkmers(x1)+...+nkmers(xm) >= KBIG,
     //     nkmers(x1)+...+nkmers(xm-1) < KBIG,
     //     nkmers(y1)+...+nkmers(yn) >= KBIG,
     //     nkmers(y1)+...+nkmers(yn-1) < KBIG.

     for ( int u = 0; u < nuni; u++ )
     {    if ( nkmers[u] < stud_size ) continue;

          // Iterate until we get somewhere.

          int KK;
          set< vec<int> > unitsd;
          digraphV< vec<int> > M;
          for ( KK = KBIG; KK >= 300; KK -= 100 )
          {    if (VERBOSE) cout << "u = " << u << ", trying KK = " << KK << endl;

               // Find initial units.

               vec< vec<int> > units0;
               for ( int j = 0; j < cindex1[u].isize( ); j++ )
               {    int id = cindex1[u][j].first, p = cindex1[u][j].second;
                    const vec<int>& y = ucores0[id];
                    int n = 0, l;
                    for ( l = p + 1; l < y.isize( ); l++ )
                    {    n += nkmers[ y[l] ];
                         if ( n >= KK ) break;    }
                    if ( l == y.isize( ) ) continue;
                    vec<int> z;
                    z.SetToSubOf( y, p, l - p + 1 );
                    if ( !Member( units0, z ) ) units0.push_back(z);    }

               // Extend units.

               vec< vec<int> > unitsp, unitsf;
               unitsd.clear( );
               for ( int j = 0; j < units0.isize( ); j++ )
               {    if ( nkmers[ units0[j].back( ) ] >= stud_size )
                         unitsf.push_back( units0[j] );
                    else unitsp.push_back( units0[j] );
                    unitsd.insert( units0[j] );    }
               while( unitsp.nonempty( ) )
               {    vec<int> x = unitsp.back( );
                    unitsp.pop_back( );
                    int d, nd = 0;
                    for ( d = 0; d < x.isize( ); d++ )
                    {    nd += nkmers[ x[d] ];
                         if ( nd >= KK ) break;    }
                    for ( int j = 0; j < cindex1[ x[d+1] ].isize( ); j++ )
                    {    int id = cindex1[ x[d+1] ][j].first;
                         int p = cindex1[ x[d+1] ][j].second;
                         const vec<int>& y = ucores0[id];
                         Bool mismatch = False;
                         for ( int l = d + 1; l < x.isize( ); l++ )
                         {    if ( l - (d+1) + p >= y.isize( ) 
                                   || x[l] != y[ l - (d+1) + p ] )
                              {    mismatch = True;
                                   break;    }    }
                         if (mismatch) continue;
                         int n = 0, r;
                         for ( r = x.isize( ) - (d+1) + p; r < y.isize( ); r++ )
                         {    n += nkmers[ y[r] ];
                              if ( n >= KK ) break;    }
                         if ( n < KK ) continue;
                         vec<int> z;
                         z.SetToSubOf( y, p, r - p + 1 );
                         if ( Member( unitsd, z ) ) continue;
                         if ( nkmers[ z.back( ) ] >= stud_size ) unitsf.push_back(z);
                         else unitsp.push_back(z);
                         unitsd.insert(z);    }    }

               if (VERBOSE)
               {    cout << "\nu = " << u << ", found " << unitsd.size( ) 
                         << " units, ";
                    PRINT( unitsf.size( ) );
                    vec<int> studs;
                    for ( int i = 0; i < unitsf.isize( ); i++ )
                         studs.push_back( unitsf[i].back( ) );
                    UniqueSort(studs);
                    cout << "unitsf reaches studs:";
                    for ( int i = 0; i < studs.isize( ); i++ )
                         cout << " " << studs[i];
                    cout << "\n";    }

               // Build a graph M from the units.  To do this we define the two 
               // halves (h1, h2) of each unit.

               vec< vec<int> > verts;
               for ( set< vec<int> >::iterator i = unitsd.begin( ); 
                    i != unitsd.end( ); i++ )
               {    verts.push_back(*i);    }
               int nv = verts.size( );
               vec< vec<int> > h1, h2;
               for ( int j = 0; j < nv; j++ )
               {    const vec<int>& x = verts[j];
                    int n = 0, l;
                    for ( l = 0; l < x.isize( ); l++ )
                    {    n += nkmers[ x[l] ];
                         if ( n >= KK ) break;    }
                    vec<int> x1, x2;
                    x1.SetToSubOf( x, 0, l+1 );
                    x2.SetToSubOf( x, l+1, x.isize( ) - (l+1) );
                    h1.push_back(x1), h2.push_back(x2);    }
               vec<int> ids1( nv, vec<int>::IDENTITY );
               SortSync( h1, ids1 );
               vec< vec<int> > from(nv), to(nv);
               for ( int v1 = 0; v1 < nv; v1++ )
               {    int low = LowerBound( h1, h2[v1] ); 
                    int high = UpperBound( h1, h2[v1] );
                    for ( int j = low; j < high; j++ )
                    {    int v2 = ids1[j];
                         from[v1].push_back(v2), to[v2].push_back(v1);    }    }
               for ( int v = 0; v < nv; v++ )
               {    Sort( from[v] ), Sort( to[v] );    }
               M.Initialize( from, to, verts );

               // Finds M's sources, which are units starting with u, and 'good'
               // sinks, which are units ending at a stud that is linked to from u.

               vec<int> sources, sinks;
               for ( int j = 0; j < verts.isize( ); j++ )
               {    if ( verts[j].front( ) == u ) sources.push_back(j);
                    int v = verts[j].back( );
                    if ( nkmers[v] >= stud_size && Member( G.From(u), v ) )
                         sinks.push_back(j);    }
     
               // Find all vertices that lie between a source and a good sink.  Then
               // delete all the other vertices.

               vec<int> s1, s2, t;
               M.GetSuccessors(sources, s1), M.GetPredecessors(sinks, s2);
               for ( int v = 0; v < nv; v++ )
                    if ( !BinMember( s1, v ) || !BinMember( s2, v ) ) t.push_back(v);
               M.DeleteVertices(t);
               if (VERBOSE) PRINT( M.N( ) );
               if ( VERBOSE && !M.Acyclic( ) ) cout << "cyclic!\n";
               
               // OK?

               if ( M.N( ) > 0 ) break;    }

          // Print the graph.

          if (VERBOSE)
          {    cout << "\n";
               // M.PrettyDOT( cout, False, True );
                    }

          // Find all paths through the graph (if it is acyclic and there are not 
          // too many paths) and delete the ones having very low coverage.

          vec< vec<int> > paths, pathsx;
          Bool have_paths = False;
          if ( M.Acyclic( ) )
          {    const int max_paths = 100;
               const int cov_ratio = 5;
               if ( M.AllPaths( -1, -1, paths, max_paths ) && paths.nonempty( ) )
               {    have_paths = True;
                    vec<int> mc( paths.size( ) );
                    for ( int pi = 0; pi < paths.isize( ); pi++ )
                    {    const vec<int>& p = paths[pi];
                         vec<int> x = M.Vert( p[0] );
                         for ( int j = 1; j < p.isize( ); j++ )
                         {    const vec<int>& y = M.Vert( p[j] );
                              int d, nd = 0;
                              for ( d = 0; d < y.isize( ); d++ )
                              {    nd += nkmers[ y[d] ];
                                   if ( nd >= KK ) break;    }
                              for ( int l = d + 1; l < y.isize( ); l++ )
                                   x.push_back( y[l] );    }
                         pathsx.push_back(x);
                         vec< vec<int> > ids( x.size( ) );
                         for ( int p1 = 0; p1 < x.isize( ); p1++ )
                         {    int p2, n = 0;
                              for ( p2 = p1; p2 < x.isize( ); p2++ )
                              {    n += nkmers[ x[p2] ];
                                   if ( n >= KK ) break;    }
                              if ( n < KK ) continue;
                              vec<int> z;
                              z.SetToSubOf( x, p1, p2 - p1 + 1 );
                              for ( int j = 0; j < cindex1[ x[p1] ].isize( ); j++ )
                              {    int id = cindex1[ x[p1] ][j].first;
                                   int q = cindex1[ x[p1] ][j].second;
                                   if ( ucores0[id].Contains( z, q ) )
                                   {    for ( int l = p1; l <= p2; l++ )
                                             ids[l].push_back(id);    }    }    }
                         vec<int> cov( x.size( ), 0 );
                         for ( int j = 0; j < x.isize( ); j++ )
                         {    UniqueSort( ids[j] );
                              cov[j] = ids[j].size( );    }
                         mc[pi] = Min(cov);
                         if (VERBOSE)
                         {    cout << "path " << pi+1 << ", mc = " << mc[pi] << ": ";
                              x.Print(cout);
                              cout << "\n";    }    }
                    int mcm = Max(mc);
                    vec<int> keep, to_delete;
                    vec<Bool> paths_to_delete( paths.size( ), False );
                    for ( int pi = 0; pi < paths.isize( ); pi++ )
                    {    if ( mc[pi] * cov_ratio > mcm ) keep.append( paths[pi] );
                         else paths_to_delete[pi] = True;    }
                    EraseIf( paths, paths_to_delete );
                    EraseIf( pathsx, paths_to_delete );
                    UniqueSort(keep);
                    for ( int j = 0; j < M.N( ); j++ )
                         if ( !BinMember( keep, j ) ) to_delete.push_back(j);
                    if (VERBOSE)
                    {    cout << "\ndeleting " << to_delete.size( ) << " vertices" 
                              << endl;    }
                    M.DeleteVertices(to_delete);    }    }
          if (VERBOSE)
          {    cout << "\n";
               for ( int i = 0; i < M.N( ); i++ )
               {    cout << "vertex " << i << " = ";
                    M.Vert(i).Print(cout);
                    cout << "\n";    }    }

          // Check for existence of paths within accepted.

          vec<int> s, sp;
          for ( int j = 0; j < aindex[u].isize( ); j++ )
               s.push_back( aindex[u][j].first );
          UniqueSort(s);
          A.GetSuccessors( s, sp );
          Bool reaches_stud = False;
          vec<int> studs;
          for ( int j = 0; j < sp.isize( ); j++ )
          {    const vec<int>& x = acc[ sp[j] ];
               int m;
               for ( m = x.isize( ) - 1; m >= 0; m-- )
                    if ( x[m] == u ) break;
               for ( int l = m + 1; l < x.isize( ); l++ )
                    if ( nkmers[ x[l] ] >= stud_size ) studs.push_back( x[l] );    }
          UniqueSort(studs);
          if (VERBOSE)
          {    cout << "\nreaches studs:";
               for ( int j = 0; j < studs.isize( ); j++ )
                    cout << " " << studs[j];
               cout << "\n";    }
          if ( studs.empty( ) )
          {    if (have_paths) 
               {    if (VERBOSE)
                         cout << "appending " << paths.size( ) << " paths" << endl;
                    exp2_verts.append(pathsx);    }
               else
               {    for ( int v = 0; v < M.N( ); v++ )
                         exp2_verts.push_back( M.Vert(v) );    }
               if (VERBOSE) cout << "\nDOES NOT REACH STUD!\n";    }

          // Find long stretches of reads that align to vertices, and define
          // coverage based on this.

          if (VERBOSE) cout << "\n";
          const int MM = KK;
          vec< vec<int> > rcov( M.N( ) );
          for ( int i = 0; i < M.N( ); i++ )
          {    const vec<int>& x = M.Vert(i);
               rcov[i].resize( x.size( ), 0 );
               vec< vec<int> > ids( x.size( ) );
               for ( int p1 = 0; p1 < x.isize( ); p1++ )
               {    int p2, n = 0;
                    for ( p2 = p1; p2 < x.isize( ); p2++ )
                    {    n += nkmers[ x[p2] ];
                         if ( n >= MM ) break;    }
                    if ( n >= MM )
                    {    vec<int> y;
                         y.SetToSubOf( x, p1, p2 - p1 + 1 );
                         for ( int j = 0; j < cindex1[ x[p1] ].isize( ); j++ )
                         {    int id = cindex1[ x[p1] ][j].first;
                              int q = cindex1[ x[p1] ][j].second;
                              if ( ucores0[id].Contains( y, q ) )
                              {    for ( int l = p1; l <= p2; l++ )
                                        ids[l].push_back(id);    }    }    }    }
               for ( int j = 0; j < x.isize( ); j++ )
               {    UniqueSort( ids[j] );
                    rcov[i][j] = ids[j].size( );    }
               if (VERBOSE)
               {    cout << "coverage of vertex " << i << " = ";
                    rcov[i].Print(cout);
                    cout << "\n";    }    }    }

     flush(cout);
     if ( !EXP2_SAVE ) Scram();
     }

     if (EXP2_SAVE) accepted.append(exp2_verts);
                    
     // Assess results.

     if (VALIDATE)
     {    cout << "\n" << nwarnings << " warnings" << endl;
          Assess( accepted, genome, genome2, unibases, K, LG, Glocs, uperfect );    }

     // Do something with results.

     vecbasevector all;
     if (WRITE)
     {    cout << "\n" << Date( ) << ": sorting results" << endl;
          UniqueSort(accepted);
          cout << Date( ) << ": building all" << endl;
          for ( int j = 0; j < accepted.isize( ); j++ )
          {    const vec<int>& x = accepted[j];
               basevector b = unibases[ x[0] ];
               for ( int l = 1; l < x.isize( ); l++ )
               {    b.resize( b.isize( ) - (K-1) );
                    b = Cat( b, unibases[ x[l] ] );    }
               all.push_back_reserve(b);    }    }
  
     // If nothing new was found then at least write out initial unibases

     if ( all.size() == 0 )
       all = unibases;

     // Write output.

     String OUT_HEAD = "extended.long";
     String outhead = run_dir + "/" + OUT_HEAD;
     if (WRITE) 
     {    cout << Date( ) << ": building output" << endl;
          vecKmerPath newpaths, newpathsrc, newunipaths;
          vec<tagged_rpint> newpathsdb, newunipathsdb;
          ReadsToPathsCoreY( all, KBIG, newpaths, newpathsrc, newpathsdb,
               run_dir + "/LRC", NUM_THREADS );
          Unipath( newpaths, newpathsrc, newpathsdb, newunipaths, newunipathsdb );
	  digraph A;
	  BuildUnipathAdjacencyGraph( newpaths, newpathsrc, newpathsdb, 
               newunipaths, newunipathsdb, A);
          KmerBaseBroker newkbb( KBIG, newpaths, newpathsrc, newpathsdb, all );
          vecbasevector newunibases;
          newunibases.reserve(newunipaths.size());
          for ( size_t i = 0; i < newunipaths.size( ); i++ )
               newunibases.push_back( newkbb.Seq( newunipaths[i] ) );
          cout << Date( ) << ": writing output files" << endl;
          String KBIGS = ToString(KBIG);
          newunipaths.WriteAll( outhead + ".unipaths.k" + KBIGS );
          BinaryWriter::writeFile( 
               outhead + ".unipathsdb.k" + KBIGS, newunipathsdb );
	  BinaryWriter::writeFile( outhead + ".unipath_adjgraph.k" + KBIGS, A );
          newunibases.WriteAll( outhead + ".unibases.k" + KBIGS );    }
     cout << Date( ) << ": done" << endl;    }
