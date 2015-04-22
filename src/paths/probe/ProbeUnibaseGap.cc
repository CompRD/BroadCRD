///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/**
 * Determing the valid unibase gaps 
*/

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "Vec.h"
#include "VecUtilities.h"
#include "math/Functions.h"
#include "paths/PdfEntry.h"
#include "paths/UnibaseUtils.h"
#include "paths/UnibaseCopyNumber3Core.h"

#include "paths/UnibaseCopyNumber3LowerCN.h"
#include "paths/LinkingPairs.h"


// A simple lass for probability function defined in range [xmin, xmax)
class IntProb {
 private:
  int min, max; 
  vec<double> p; // probability distribution
 public:
  IntProb ( const vec<double>& array, const int MIN_SEP = 0 ) { InitFromArray( array, MIN_SEP); };
  IntProb(){};

  void InitFromArray( const vec<double>& array, const int MIN_SEP = 0 ) {
       int len = array.size();
       ForceAssertGt( len, 0 );
       p = array;
       min = MIN_SEP;
       max = min + len;
  };

  // the implementation of virtual methods
  double p0(int t) const { return seekArray(t, p); }
  double L() const  { return min;  }
  double U() const  { return max;  }
 private:
  double seekArray(int t, const vec<double>& dd) const {
    int index = t - min;
    if ( index < 0 || index >= dd.isize() ) return 0;
    return dd[index];
  }
};


// obtain the len-bias corrected distribution
void GetIntProbs ( 
          const LinkingPairs& linking,
          vec< IntProb >& int_distrs
          ) 
{
     int nlibs = linking.NLibs();
     int_distrs.assign( nlibs, IntProb() );
     for ( int libid = 0; libid < nlibs; libid++ ) {
          const StdMap<int,int>& dist = linking.GetDists(libid);
          //cout << "histx" << libid << " ";
          //for( map<int,int>::const_iterator it = dist.begin(); it != dist.end(); it++ )
          //     cout << it->first << " ";
          //cout << endl;
          //cout << "histy" << libid << " ";
          //for( map<int,int>::const_iterator it = dist.begin(); it != dist.end(); it++ )
          //     cout << it->second << " ";
          //cout << endl;
          const vec<bool>& unibase_pick = linking.GetTigPick( libid );
          vec<int> lens;
          for ( int i = 0; i < linking.NTigs(); i++ )
               if ( unibase_pick[i] ) lens.push_back( linking.GetLen(i) );
          Sort(lens);

          vec<double> distr_array;
          int cutoff_lb ;
          GetDistribution( dist, lens, distr_array, cutoff_lb,  true);
          //cout << scientific; 
          //for ( int i = 0; i < distr_array.isize(); i++ ) {
          //     if ( i % 10 != 0 ) continue;
          //     int x = i + cutoff_lb;
          //     cout << "distr" << libid << " " << x << " " << distr_array[i] << endl;
          //}
          int_distrs[libid].InitFromArray( distr_array, cutoff_lb );
	  //distrs[libid].FromArray( distr_array, cutoff_lb );
          //for ( int x = int_distrs[libid].L(); x < int_distrs[libid].U(); x++ ) {
          //     if ( x % 10 != 0 ) continue;
          //     cout << "distr" << libid << " " << x << " " << int_distrs[libid].p0(x) << endl;
          //}
     }
}

double UnibaseGap3( int gap, const vec< IntProb >& distrs, 
          int u1, int u2,
          int len1, int len2,
         const vec< vec< pair<int,int> > >& links12,
         const LibCov_t & starts1,
         const LibCov_t & stops2,
         // end of input
	 double & prob_pos,
	 double & prob_neg
    )
{
     double linking_rate = 0.99;   // Chance forming a link
     double noise_rate = 1e-10;    // the chance of forming a random link

     int nlibs = distrs.size();
     ForceAssertEq( nlibs, links12.isize() );
     ForceAssertEq( nlibs, starts1.isize() );
     ForceAssertEq( nlibs, stops2.isize() );
     prob_pos = 0;
     prob_neg = 0;
     for ( int ilib = 0; ilib < nlibs; ilib++ ) {
          int max = distrs[ilib].U();
          multiset< LinkEnd_t >::iterator it = starts1[ilib].begin();
          for ( ; it !=  starts1[ilib].end(); it++ ) {
               int start_pos = it->first;
               int range1 = 1;
               int end_tig = it->second;
               int end_pos = it->third;
               int link_end_pos = -1;
               if ( end_tig == u2 ) link_end_pos = end_pos;  // mark the end position
               // the other possible links
               multiset< LinkEnd_t >::iterator it2 = stops2[ilib].begin();
               for ( ; it2 !=  stops2[ilib].end(); it2++ ) {
                    int end_pos = it2->first;
                    int range2 = 1;
                    int dist = len1 - start_pos + gap + end_pos;
                    double p_dist = distrs[ilib].p0(dist) * range2;
                    if ( end_pos != link_end_pos )
                         prob_neg += log( 1 - p_dist * linking_rate - noise_rate );
                    else
                         prob_pos += log( p_dist + noise_rate ); 
               }
          }
     }
     double score = prob_pos + prob_neg;
     return score;
}


double UnibaseGap4( int gap, const vec< IntProb >& distrs, 
	            int u1, int u2,
	            int len1, int len2,
                    const vec< vec< pair<int,int> > >& links12,
                    const vec< StdMap<int,double> > & starts1,
                    const vec< StdMap<int,double> > & stops2,
                    // end of input
	            double & prob_pos,
	            double & prob_neg 
    )
{
     //const double linking_rate = 0.75;   // Chance forming a link
     //const double linking_rate = 0.62;   // Chance forming a link
     const double linking_rate = 0.97;   // Chance forming a link
     const double noise_rate = 1e-10;    // the chance of forming a random link

     int nlibs = distrs.size();
     ForceAssertEq( nlibs, links12.isize() );
     ForceAssertEq( nlibs, stops2.isize() );
     prob_pos = 0;
     prob_neg = 0;
     for ( int ilib = 0; ilib < nlibs; ilib++ ) {
          int dist_lb = distrs[ilib].L();
          int dist_ub = distrs[ilib].U();
          map<int, double>::const_iterator it = starts1[ilib].begin();
          for ( ; it !=  starts1[ilib].end(); it++ ) {
               int start_pos = it->first;
	       double cov1 = it->second;
	       //if ( len1 - start_pos + gap + len2 < dist_lb ) continue;
	       //if ( len1 - start_pos + gap + 0 > dist_ub ) break;
               map< int,double>::const_iterator it2 = stops2[ilib].begin();
               for ( ; it2 !=  stops2[ilib].end(); it2++ ) {
                    int end_pos = it2->first;
                    int dist = len1 - start_pos + gap + end_pos;
		    if ( dist < dist_lb ) continue;
		    if ( dist > dist_ub ) break;
		    double cov2 = it2->second;
		    //double cov = ( cov1 + cov2 ) * 0.5;
		    //double cov = sqrt( cov1 * cov2 );
		    double cov = Min( cov1 , cov2 );
                    double p_dist = Min( 1.0, ( noise_rate + distrs[ilib].p0(dist) )* cov * linking_rate );
		    prob_neg += log( 1 - p_dist );
               }
          }
	  for( int k = 0; k < links12[ilib].isize(); k++ ) {
	       int start_pos = links12[ilib][k].first;
	       int end_pos = links12[ilib][k].second;
	       int dist = len1 - start_pos + gap + end_pos;
	       double cov1 = 0, cov2 = 0;
	       map<int, double>::const_iterator it = starts1[ilib].find( start_pos );
	       if ( it != starts1[ilib].end() ) cov1 = it->second;
	       it =  stops2[ilib].find( end_pos );
	       if ( it != stops2[ilib].end() ) cov2 = it->second;
	       //double cov = ( cov1 + cov2 ) * 0.5 ;
	       //double cov = sqrt( cov1 * cov2 );
	       double cov = Min( cov1 , cov2 );
	       //double p_dist = Min( 1.0, distrs[ilib].p0(dist) * cov * linking_rate );
	       //if ( p_dist < noise_rate ) p_dist = noise_rate;
	       double p_dist = Min( 1.0, ( noise_rate + distrs[ilib].p0(dist)) * cov * linking_rate );
	       prob_pos += log( p_dist ); 
	       prob_neg -= log( 1 - p_dist );
	  }
     }
     double score = prob_pos + prob_neg;
     return score;
}


// Experimental code that find all links between unibases
// and determine the most probable gap size and standard deviation between them
void GenerateUnibaseLinks( 
    // input
    const String & run_dir,
    const String & JUMP_READS,
    const int K,
    const vec<int> & to_rc,
    const vec<segalign> & SUGS,
    const vecbasevector & unibases,
    const vec<double>&  CN_raw,
    // output
    LinkingPairs& linking,
    // parameters
    bool VERBOSITY = 0
    )
{
     String head = run_dir + "/" + JUMP_READS;
     PairsManager pairs( head + ".pairs" );
     int nlibs = pairs.nLibraries();
     int ntigs = unibases.size();

     // For each pair of unipaths that are connected by two or more links, predict
     // their order and separation.  Both orders may be possible.  Note that a more 
     // wholistic approach may be needed.
     cout << Date( ) << ": set up for link computation" << endl;
     vecbasevector reads( head + ".fastb" );
     uint64_t nreads = reads.size( );
     vec<size_t> S_START_RID(nreads+1);
     {    size_t SEG_POS = 0;
          for ( uint64_t rid = 0; rid <= nreads; rid++ )
          {   while( SEG_POS < SUGS.size( ) && SUGS[SEG_POS].rid < rid ) ++SEG_POS;
              S_START_RID[rid] = SEG_POS;    }    }
     cout << Date( ) << ": compute links from "
         << ToStringAddCommas( pairs.nPairs( ) ) << " pairs" << endl;
     
     // prepare the alignments for each read. must be uniquely aligned
     vec< pair<int,int> > aligns( nreads, make_pair(-1,-1) );
     for ( size_t i = 0; i < pairs.nPairs( ); i++ ) {    
          longlong id1 = pairs.ID1(i), id2 = pairs.ID2(i);
          int libid = pairs.libraryID(i);
          // Only when reads are uniquely alignment !
          if ( S_START_RID[id1+1] - S_START_RID[id1] == 1 ) {
               size_t j1 = S_START_RID[id1];
               int x1 = SUGS[j1].upos - SUGS[j1].rpos; 
               int u1 = SUGS[j1].u;
	       int nu1 = unibases[u1].size( );
	       if ( x1 >= 0 && x1 < nu1 )
                         aligns[ id1 ] = make_pair( u1, x1 );
          }
          if ( S_START_RID[id2+1] - S_START_RID[id2] == 1 ) {
               size_t j2 = S_START_RID[id2];
               int x2 = SUGS[j2].upos - SUGS[j2].rpos; 
               int u2 = SUGS[j2].u;
	       int nu2 = unibases[u2].size( );
	       if ( x2 >= 0 && x2 < nu2 )
                    aligns[ id2 ] = make_pair( u2, x2 );
          }
     }
     // Gather the linking infomation
     vec< int > unibase_lens( unibases.size(), -1); 
     for ( size_t u = 0; u < unibases.size(); u++ )
          unibase_lens[u] = unibases[u].size();
     cout << "aligns.size()= " << aligns.size() << endl;

     linking.Init( nlibs, unibase_lens);
     GatherLinks( pairs, aligns, to_rc, linking );

     if ( VERBOSITY < 3 ) return;

     // the distribution
     vec< IntProb > int_distrs( nlibs);
     GetIntProbs ( linking, int_distrs );

     // check the model
     #pragma omp parallel for
     for( int itig = 0; itig < ntigs; itig += 1 ) {
          int len = unibases[itig].size();
          const LibCov_t & starts = linking.GetStarts( itig );
          const LibCov_t & stops = linking.GetStops( itig );
          if ( starts[0].empty() || stops[0].empty() ) continue;
          if ( len < 5000 ) continue;
	  // find the coverage
	  vec< StdMap<int,double> > starts1(nlibs);
	  vec< StdMap<int,double> > stops2(nlibs);
	  linking.GetStartsSampled(itig, starts1);
	  linking.GetStopsSampled(itig, stops2);
	  //int delta = 50;
	  //linking.GetStartsSmoothed(itig, starts1, delta);
	  //linking.GetStopsSmoothed(itig, stops2, delta);
	  // what is the integration of the model over the contig?
	  // cov(x) * cov(y) * p( y - x)
	  int libid = 0;
	  double npred1 = 0;
	  double npred2 = 0;
	  double npred3 = 0;
	  for ( auto it1 = starts1[libid].begin();
		    it1 != starts1[libid].end(); it1++ )
	  for ( auto it2 = stops2[libid].begin();
		    it2!= stops2[libid].end(); it2++ )
	  {
	       int dist = it2->first - it1->first ;
	       if ( dist < 4000 ) continue;
	       if ( dist > 5000 ) break;
	       //double s = 1;
	       double s = it1->second * it2->second;
	       double f = int_distrs[libid].p0( dist ) / s;
	       npred1 += f *  sqrt(it1->second * it2->second) ;
	       npred2 += f * (it1->second + it2->second) * 0.5 ;
	       npred3 += f * Min( it1->second , it2->second ) ;
	  }
          int local_links = 0;
          int other = 0;
          int miss = 0;
	  int sel_links = 0;
          for ( multiset< LinkEnd_t >:: iterator it = starts[0].begin(); it != starts[0].end(); it++ ) {
               if ( it->second == itig ) local_links++;
               if ( it->second == -1 ) miss++  ;
               if ( it->second != itig && it->second != -1 ) other++;
               if ( it->second == itig && it->third - it->first < 5000
			               && it->third - it->first > 4000 ) sel_links++;
          }
          #pragma omp critical
	  {
               cout << "itig= " << itig << " ";
               cout << "copy= " << CN_raw[itig] << " ";
               cout << "len= " << len << " ";
               cout << "nstarts= " << starts[0].size() << " ";
               cout << "nstops= " << stops[0].size() << " ";
               cout << "npred1= " << npred1 << " ";
               cout << "npred2= " << npred2 << " ";
               cout << "npred3= " << npred3 << " ";
               cout << "sel_links= " << sel_links << " ";
               cout << "nlocal= " << local_links << " ";
               cout << "nmiss= " << miss << " ";
               cout << "nother= " << other << endl;
          }
     }
}


void EstimateUnibaseGap(
    // input
    const int U1,
    const int U2,
    const LinkingPairs& linking, 
    // output
    vec<ulink_with_uids> & condensed_links, 
    // parameters
    bool VERBOSITY = 0
    )
{
     cout << Date( ) << ": condense links" << endl;
     int nlibs = linking.NLibs();
     int ntigs = linking.NTigs();
     const int min_links_initial = 2;
     int dev_mult = 3;

     // the distribution
     vec< IntProb > int_distrs( nlibs);
     GetIntProbs ( linking, int_distrs );

     cout << scientific;
     cout << setprecision(3);

     //// print distribution
     //for ( int libid = 0; libid < nlibs; libid++ ) {
     //     IntProb& pb = int_distrs[libid];
     //     for ( int x = pb.L(); x < pb.U(); x++ )
     //          cout << "lib" << libid << " " << x << " " <<  pb.p0(x) << endl;
     //}
     //return;

     vec< pair<int,int> > all_linked = linking.GetAllLinked( );
     for ( int i = 0; i < all_linked.isize(); i++ ) {
          // We get all the linked pairs where there are read pairs that start from (rc align)
          // u1 and stop at (fw align) u2. Two possible arrangements :
          // (1)
          //                  -----------------u1-------------- (gap) ------u2------
          //                                         <---@ r1             @-->r2     (jump)
          // (2)
          // ----u2---- (gap) -----------------u1--------------                      
          //  @-->r2              <---@ r1                                           (non-jump)
          // Case (2) can be treated as case (1) with a large negative gap.
          //
          //
          // We should also check pairs which start from u2 and end at u1.
          // (1)
          //                  -----------------u1-------------- (gap) ------u2------ 
          //                                         @---> r1             <--@r2    (non-jump)
          // (2)
          // ----u2---- (gap) -----------------u1--------------
          //  <--@r2              @---> r1                                          (jump)

          int u1 = all_linked[i].first;
          int u2 = all_linked[i].second;
          if ( u1 == u2 ) continue;
          const LibLink_t & ulinks12= linking.GetLinks( u1, u2 ); // starts from u1, ends at u2
          const LibLink_t & ulinks21= linking.GetLinks( u2, u1 ); // starts from u2, ends at u1
          int nlink_jp = 0, nlink_nj = 0;  // jump and nonjump links
          for ( int libid = 0; libid < nlibs; libid++ ) {
               nlink_jp += ulinks12[libid].size();
               nlink_nj +=  ulinks21[libid].size();
          }
          int nlinks = nlink_jp + nlink_nj; 
          //if ( nlinks < 20 ) continue;

          if ( u1 != U1|| u2 != U2) continue;
          int len1 = linking.GetLen(u1);
          int len2 = linking.GetLen(u2);

	  // find the coverage
	  vec< StdMap<int,double> > starts1(nlibs);
	  vec< StdMap<int,double> > starts2(nlibs);
	  vec< StdMap<int,double> > stops1(nlibs);
	  vec< StdMap<int,double> > stops2(nlibs);
	  int delta = 50;
	  linking.GetStartsSmoothed(u1, starts1, delta);
	  linking.GetStartsSmoothed(u2, starts2, delta);
	  linking.GetStopsSmoothed(u1, stops1, delta);
	  linking.GetStopsSmoothed(u2, stops2, delta);
	  //linking.GetStartsSampled(u1, starts1 );
	  //linking.GetStartsSampled(u2, starts2);
	  //linking.GetStopsSampled(u1, stops1);
	  //linking.GetStopsSampled(u2, stops2);
	  DumpLinkingInfo( cout, starts1, stops2, ulinks12, 10, true, len1 );

          vec< pair<double, int> > results;
	  vec< vec<double> > scores_comp;
          #pragma omp parallel for
          for ( int sep = -2000; sep <= 8000; sep += 100 ) {
               //double prob = UnibaseGap3 ( sep, int_distrs, 
               //               u1, u2,
               //               len1,len2,
               //               ulinks12, linking.GetStarts(u1), linking.GetStops(u2), 
               //               ulinks21, linking.GetStarts(u2), linking.GetStops(u1) );
	       vec<double> scores(4,0);
               double prob = UnibaseGap4 ( sep, int_distrs, 
                              u1, u2,
                              len1,len2,
                              ulinks12, starts1, stops2, scores[0], scores[1] ) ;
	       prob +=       UnibaseGap4 ( - sep - len1 - len2, int_distrs, 
                              u2, u1,
                              len2,len1,
                              ulinks21, starts2, stops1, scores[2], scores[3] ) ;
               #pragma omp critical
               { results.push( prob, sep ); 
                 scores_comp.push_back( scores );
	       }
          }

          #pragma omp critical
          {
          ReverseSortSync( results, scores_comp );
          double max_prob = results.front().first;
          double max_prob_gap = results.front().second;
          if ( VERBOSITY >= 1 ) 
	  {
               for( int i = 0; i < results.isize(); i++ )
                    cout << "prob " << results[i].second << " " << pow(10, results[i].first - max_prob) 
                         << " " << results[i].first
			 << " " << scores_comp[i][0] 
			 << " " << scores_comp[i][1] 
			 << " " << scores_comp[i][2] 
			 << " " << scores_comp[i][3] 
			 <<  endl;

               for ( int libid = 0; libid < nlibs; libid++ ) {
                    for ( int k = 0; k < ulinks12[libid].isize(); k++ ) 
                         cout << "link12 " << libid << " " << ulinks12[libid][k].first << " "<< ulinks12[libid][k].second << endl;
               }
               for ( int libid = 0; libid < nlibs; libid++ ) {
                    for ( int k = 0; k < ulinks21[libid].isize(); k++ ) 
                         cout << "link21 " << libid << " " << ulinks21[libid][k].first << " " << ulinks21[libid][k].second << endl;
               }
               cout << "pair " << u1 << "_" << u2 << " " << "nlinks = " << nlinks << 
                    " ( " << nlink_jp << " u1->u2, " << nlink_nj << " u2->u1 ) " << endl;
	       cout << "len1= " << len1 << " nstarts1= " ;
               for ( int libid = 0; libid < nlibs; libid++ ) cout << linking.GetNStarts(libid, u1) << " ";
	       cout << endl;
	       cout << "len2= " << len2 << " nstops2= " ;
               for ( int libid = 0; libid < nlibs; libid++ ) cout << linking.GetNStops(libid, u2) << " ";
	       cout << endl;
	        
               cout << "estimated_gap= " << max_prob_gap << endl;
	  } 
          //if ( sep + max_devs * dev <  -(K-1) ) continue;
          //int start1 = 0, stop2 = 0; // those numbers are actually not used afterwards
          //condensed_links.push( u1, u2, sep, dev, start1, stop2, nlinks );
          }
     }
}


void DisplayUnibaseGap(
    // input
    const int U1,
    const int U2,
    const LinkingPairs& linking, 
    // parameters
    bool VERBOSITY = 0
    )
{
     int u1 = U1, u2 = U2;
     int nlibs = linking.NLibs();
     int ntigs = linking.NTigs();
     int len1 = linking.GetLen(u1);
     int len2 = linking.GetLen(u2);
     const LibLink_t & ulinks12= linking.GetLinks( u1, u2 ); // starts from u1, ends at u2
     const LibLink_t & ulinks21= linking.GetLinks( u2, u1 ); // starts from u2, ends at u1
     // find the coverage
     vec< StdMap<int,double> > starts1(nlibs);
     vec< StdMap<int,double> > starts2(nlibs);
     vec< StdMap<int,double> > stops1(nlibs);
     vec< StdMap<int,double> > stops2(nlibs);
     int delta = 50;
     linking.GetStartsSmoothed(u1, starts1, delta);
     linking.GetStartsSmoothed(u2, starts2, delta);
     linking.GetStopsSmoothed(u1, stops1, delta);
     linking.GetStopsSmoothed(u2, stops2, delta);
     DumpLinkingInfo( cout, starts1, stops2, ulinks12, 10, true, len1 );
     const LibCov_t & ss1 = linking.GetStarts( u1 );
     const LibCov_t & ss2 = linking.GetStarts( u2 );

     if ( VERBOSITY >= 1 ) 
     {
          int nlink_jp = 0, nlink_nj = 0, nlinks = 0;
          for ( int libid = 0; libid < nlibs; libid++ ) {
               //for ( int k = 0; k < ulinks12[libid].isize(); k++ ) 
               //     cout << "link12 " << libid << " " << ulinks12[libid][k].first << " "<< ulinks12[libid][k].second << endl;
               for ( multiset< triple<int,int,int> >::const_iterator it = ss1[libid].begin();
                  it != ss1[libid].end(); it++ ) {
                    if ( it->second == u2 )
                         cout << "link12 " << libid << " " << it->first << " "<< it->third << endl;
               }
          }
          for ( int libid = 0; libid < nlibs; libid++ ) {
               //for ( int k = 0; k < ulinks21[libid].isize(); k++ ) 
               //     cout << "link21 " << libid << " " << ulinks21[libid][k].first << " " << ulinks21[libid][k].second << endl;
               for ( multiset< triple<int,int,int> >::const_iterator it = ss2[libid].begin();
                  it != ss2[libid].end(); it++ ) {
                    if ( it->second == u1 )
                         cout << "link21 " << libid << " " << it->first << " "<< it->third << endl;
               }
          }
          cout << "pair " << u1 << "_" << u2 << " " << "nlinks = " << nlinks << 
               " ( " << nlink_jp << " u1->u2, " << nlink_nj << " u2->u1 ) " << endl;
          cout << "len1= " << len1 << " nstarts1= " ;
          for ( int libid = 0; libid < nlibs; libid++ ) cout << linking.GetNStarts(libid, u1) << " ";
          cout << endl;
          cout << "len2= " << len2 << " nstops2= " ;
          for ( int libid = 0; libid < nlibs; libid++ ) cout << linking.GetNStops(libid, u2) << " ";
          cout << endl;
     } 
}
int main( int argc, char** argv ) {
  
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_Int(K);
  
  // Infixes for input/output file names.
  CommandArgument_String( READS );
  CommandArgument_String_OrDefault( UNIBASES, "unibases" );
  
  // Runtime control.

  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");

  CommandArgument_String_OrDefault(JUMP_READS, "jump_reads_filt_cpd");
  CommandArgument_Int_OrDefault(MAX_PLACEMENTS, 50);
  CommandArgument_Int_OrDefault(VERBOSITY, 0);

  CommandArgument_Int(U1);
  CommandArgument_Int(U2);
  CommandArgument_Bool_OrDefault(ESTIMATE_GAP, False);
  
  EndCommandArguments;

  // Thread control
   
  NUM_THREADS = configNumThreads(NUM_THREADS);
  omp_set_num_threads( NUM_THREADS );

  // Define directories.
  
  String data_dir = PRE + "/" + DATA;
  String run_dir = data_dir + "/" + RUN;

  String merged_head = run_dir + "/UnibaseCopyNumber3";

  // Load unibases.
  String unibases_head = run_dir + "/" + READS + "." + UNIBASES;
  
  cout << Date() << ": Loading unibases" << endl;
  String unifile = unibases_head + ".k" + ToString(K);
  vecbasevector unibases(unifile);
  cout << Date() << " End loading " << unibases.size() << " unibases" << endl;

  // Now load unipath placements on reference if needed.

  VecPlacementVec placements;
  if ( VERBOSITY >= 3 )
  {    placements.ReadAll( run_dir + "/" + READS + ".unipaths.k" + ToString(K)
      + ".locs" );    
  }

  // Load innie stats.

  vec<int> innie_sep, innie_dev;
  vec<double> innie_percent;
  fast_ifstream iin( run_dir + "/" + JUMP_READS + ".outies" );
  String line;
  while(1)
  {    getline( iin, line );
       if ( iin.fail( ) ) break;
       istrstream iline( line.c_str( ) );
       int sep, dev;
       double iper;
       iline >> sep >> dev >> iper;
       innie_sep.push_back(sep), innie_dev.push_back(dev);
       innie_percent.push_back(iper);
       if ( VERBOSITY >= 1 )
       {    cout << Date( ) << ": jump innies " << sep << " +/- " << dev << " (" 
                 << setiosflags(ios::fixed) << setprecision(1)
                 << iper << "%)" << endl;    
       }    
  }
  
  // Align jumping reads.
  
  cout << Date( ) << ": " 
       << "-----------------------------------------------------" << endl;

  // Set up ancillary data structures for unibases.
  
  size_t nuni = unibases.size( );
  vec<int> to_rc;
  UnibaseInvolution( unibases, to_rc );

  // Align reads to unibases

  cout << Date( ) << ": loading jumps" << endl;
  String temp_file = run_dir + "/tmp/temp_jumps.fastb";
  {    vecbasevector jreads( run_dir + "/" + JUMP_READS + ".fastb" );
       for ( size_t i = 0; i < jreads.size( ); i++ )
	    jreads[i].ReverseComplement( );
       cout << Date( ) << ": total jump reads = " 
	    << ToStringAddCommas( jreads.size( ) ) << endl;
       for ( size_t i = 0; i < jreads.size( ); i++ )
	    jreads[i].resize(20);
       jreads.WriteAll(temp_file);    
  }
  vec< triple<int64_t,int64_t,int> > JALIGNS;
  {    vec< triple<int64_t,int64_t,int> > jaligns;
       SearchFastb2( temp_file, unifile, 20, &jaligns, 0, MAX_PLACEMENTS );
       for ( size_t i = 0; i < jaligns.size( ); i++ )
	    if ( jaligns[i].third >= 0 ) JALIGNS.push_back( jaligns[i] );    
  }
  // Directly convert alignments into segments, then get rid of alignments.
  vec<segalign> JSEGS;
  JSEGS.resize( JALIGNS.size( ) );
  for ( size_t i = 0; i < JALIGNS.size( ); i++ ) 
  {    JSEGS[i] = segalign( True, JALIGNS[i].first, 0, JALIGNS[i].second,
	    JALIGNS[i].third );    
  }
  Destroy(JALIGNS);
  ParallelSort(JSEGS, cmp_ridx);


  // Gather and organize the linking between unibases

  cout << Date( ) << ": Generate Unibases links " << endl;
  vec<double> CN_raw( unibases.size() );
  LinkingPairs linking;
  GenerateUnibaseLinks( run_dir, JUMP_READS, K, to_rc, JSEGS, unibases, CN_raw,
            linking, VERBOSITY);

  DisplayUnibaseGap ( U1, U2, linking, VERBOSITY );

  // Esitmate gap sizes

  cout << Date( ) << ": Gap size computation" << endl;
  vec<ulink_with_uids> condensed_links;
  if ( ESTIMATE_GAP )
       EstimateUnibaseGap ( U1, U2, linking, condensed_links, VERBOSITY );

  // Compare with actual gap size

  int gap_estimate = 0;
  if ( VERBOSITY >= 2 )
  {    
       int u1 = U1;
       int u2 = U2;
       if ( placements[u1].empty( ) || placements[u2].empty( ) )
            cout << "Gap not both mapped\n";
       else
       {
            cout << "Unibase " << u1 << " are aligned to: " << endl;
            for( size_t j = 0; j < placements[u1].size(); j++ )
                 cout << placements[u1][j].Pos() << "@" <<  placements[u1][j].GenomeId()<< endl;
            cout << "Unibase " << u2 << " are aligned to: " << endl;
            for( size_t j = 0; j < placements[u2].size(); j++ )
                 cout << placements[u2][j].Pos() << "@" <<  placements[u2][j].GenomeId()<< endl;

            double infinity = 1000000000.0;
            int best_sep = 0;
            double best_offby = infinity;
            for ( size_t j1 = 0; j1 < placements[u1].size(); j1++ )
            {    for ( size_t j2 = 0; j2 < placements[u2].size(); j2++ )
     	    {    const placement& p1 = placements[u1][j1];
     		 const placement& p2 = placements[u2][j2];
     		 if ( p1.Orient( ) != p2.Orient( ) ) continue;
     		 if ( p1.GenomeId( ) != p2.GenomeId( ) ) continue;
     		 int sep;
     		 double offby;
     		 if ( p1.Fw( ) ) sep = p2.pos( ) - p1.Pos( );
     		 else sep = p1.pos( ) - p2.Pos( );
     		 offby = Abs( sep - gap_estimate );
     		 if ( offby < best_offby )
     		 {    best_sep = sep;
     		      best_offby = offby;    }    
     	    }    
            }
            cout << "Actual gap size " << best_sep << endl;
       }
  }
  return 0;    
}
