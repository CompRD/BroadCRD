///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// FixLocal.
//
// Needs AlignReads to have been run, and also SamplePairedReadDistributions.
//
// This is a very preliminary version of a program to edit an efasta assembly
// using read alignments to it.
//
// Scan assembly to locate weak spots, then locally reassemble.  The following
// steps are followed:
// 1. Look for positions providing a signal that something might be wrong.
// 2. Define windows around signal bases.
// 3. Identify reads that touch the windows.
// 4. Assemble the reads, yielding a graph.
// 5. Clean the graph.
// 6. Find all paths through the graph.
// 7. Score each path.
// 8. Pick winning paths.
// 9. Form into efasta.
//
// Known problems.
//
// 1. Presence of spurious edges.  In cases of tandem repeats (and perhaps in other
// cases), there is a tendency for the assembly graph to contain edges that are
// present only because of uncorrected sequencing errors.  There are various ways
// in which these errors might be corrected:
// (a) Running PreCorrect might help.
// (b) For each kmer x, we could test each position p on it by forming the four
// mutated kmers x_p:b, b in {A,C,G,T}.  For each b we would find all its instances
// in the corrected reads and form the list q_p:b of the associated quality scores
// at p.  Based on these data, under appropriate circumstances we might decide that
// some instances x_p:b are incorrect, and delete them from the graph.
// (c) FindErrors might be modified to do some version of (b).
//
// 2. Explosion in mapping reads.  In cases where the assembly graph contains
// enough alternatives (for example two loops may be enough), and a read contains
// many very low quality bases, mapping of the read to the graph may explode,
// because there are may be many equally likely paths for the read through the 
// graph.
//
// 3. Difficulty in scoring mappings to tandem repeats.  When a read is mapped to a
// tandem repeat, a fundamental question is the scoring of mappings that either
// terminate in a repeat copy, or leave the repeat.  Because the ends of reads tend
// to have low copy, this can be very difficult.
//
// 4. Overall computational performance.  This code has only been tested on 
// bacterial genomes.  On larger genomes the computational requirements would be
// much larger, both because the genomes are larger and because the rate of signal
// occurrence would generally be much higher.
//
// 5. Inaccurate separation estimates.  These are presumably arising because of 
// nonalignment to repeat regions.  Note also that using fragment pairs in the
// separation estimates can result in getting completely the wrong answer.
//
// 6. Not checking for signal from pairing deviation.  This would be easy.
//
// 7. Need to sort alternatives by decreasing likelihood.
//
// 8. As implemented, the bulk of run time may be spent in FindErrors.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "Qualvector.h"
#include "feudal/QualNibbleVec.h"
#include "Superb.h"
#include "efasta/EfastaTools.h"
#include "graph/Digraph.h"
#include "kmers/naif_kmer/Kmers.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "paths/AssemblyCleanupTools.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"
#include "paths/PairDistCorrection.h"
#include "paths/ReadLoc.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/Unipath.h"
#include "paths/FindErrorsCore.h"

void AddToPileup( const read_loc& rl, const basevector& b, const qualvector& q,
     const basevector& tig, vec<dumbcall>& calls )
{    align a;
     rl.GetAlign( a, b, tig );
     int p1 = a.pos1( ), p2 = a.pos2( );
     for ( int j = 0; j < a.Nblocks( ); j++ )
     {    if ( a.Gaps(j) > 0 )
          {    for ( int u = 0; u < a.Gaps(j); u++ )
                    calls[p2+u].base[4] += q[p1];
               p2 += a.Gaps(j);    }
          if ( a.Gaps(j) < 0 )
          {    for ( int u = 0; u < -a.Gaps(j); u++ )
                    calls[p2].base[5] += q[p1];
               p1 -= a.Gaps(j);    }
          for ( int x = 0; x < a.Lengths(j); x++ )
          {    calls[p2].base[ b[p1] ] += q[p1];
               ++p1; ++p2;    }    }    }

// AllSubpathTargets.  Find all paths E1-...-Em in a given HyperBasevector hb 
// that would be targets for gap-free alignment of a read of length N.  The 
// requirements for this are that
//     N <= len(E1-...-Em) and
//     N > len(E2-...-Em-1) + 1.
// Each returned path is given as a list of edge ids.

Bool AllSubpathTargets( const HyperBasevector& hb, const int N,
     vec< vec<int> >& subpaths )
{
     subpaths.clear( );
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     int K = hb.K( );
     vec< vec<int> > partials;
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
     {    vec<int> p;
          p.push_back(e);
          partials.push_back(p);    }
     while( partials.nonempty( ) )
     {    vec<int> p = partials.back( );
          partials.pop_back( );
          int len1 = hb.EdgeLengthBases( p[0] );
          for ( int j = 1; j < p.isize( ); j++ )
               len1 += hb.EdgeLengthBases( p[j] ) - (K-1);
          Bool constraint1 = ( N <= len1 );
          Bool constraint2 = True;
          if ( p.size( ) >= 2 )
          {    int len2 = hb.EdgeLengthBases( p[1] );
               for ( int j = 2; j < p.isize( ) - 1; j++ )
                    len2 += hb.EdgeLengthBases( p[j] ) - (K-1);
               if ( !( N > len2 + 1 ) ) constraint2 = False;    }
          if ( !constraint2 ) continue;
          if (constraint1) subpaths.push_back(p);
          if ( subpaths.size( ) > 1000 ) return False;
          int w = to_right[ p.back( ) ];
          for ( int j = 0; j < hb.From(w).isize( ); j++ )
          {    vec<int> q = p;
               q.push_back( hb.EdgeObjectIndexByIndexFrom( w, j ) );
               partials.push_back(q);    }    }
     return True;    }

void GapStatsAlt( vec<int> gap, vec<int> gapdev, int& gap_ave, int& gapdev_ave ){
     // If there are less than six gaps, we directly compute their mean.
     // Otherwise, we attempt to remove outliers, as follows.  We sort 
     // the gaps and extract the middle half.  From this middle half, we compute the 
     // mean and standard deviation.  Then we select those gaps lying withing 5 
     // standard deviations of the mean, and form their mean.  

     vec<NormalDistribution> S;
     if ( gap.size( ) >= 6 ){    
       vec<int> mid_gaps;
       SortSync( gap, gapdev );
       for ( unsigned int i = gap.size( )/4; i < 3*(1+gap.size( ))/4; i++ )
	 mid_gaps.push_back( gap[i] );
       float sum1 = 0, sum2 = 0;
       for ( unsigned int i = 0; i < mid_gaps.size( ); i++ ){    
	 sum1 += mid_gaps[i];
	 sum2 += float(mid_gaps[i]) * float(mid_gaps[i]);    
       }
       float n = mid_gaps.size( );
       float mean = sum1/n;
       float sd = sqrt(sum2/n - mean * mean);
       float start = mean - 5 * sd, stop = mean + 5 * sd;
       for ( unsigned int i = 0; i < gap.size( ); i++ ){    
	 if ( start <= gap[i] && gap[i] <= stop )
	   S.push( gap[i], gapdev[i] );    
       }    
     }
     else{    
       for ( int l = 0; l < gap.isize( ); l++ )
	 S.push( gap[l], gapdev[l] );    
     }
     NormalDistribution s = CombineNormalDistributions(S);
     gap_ave = int(round(s.mu_));
     gapdev_ave = int(round(s.sigma_));    
}

// Delete initial edges that do not contain the left flank, and terminal edges 
// that do no contain the right flank.

void DeleteIllegalTerminators( HyperKmerPath& h, const KmerBaseBroker& kbb,
     const basevector& left, const basevector& right )
{    while(1)
     {    vec<int> to_delete;
          for ( int v = 0; v < h.N( ); v++ )
          {    if ( h.Source(v) )
               {    for ( int j = 0; j < h.From(v).isize( ); j++ )
                    {    basevector b = kbb.Seq( h.EdgeObjectByIndexFrom( v, j ) );
                         if ( !b.ToString( ).Contains( left.ToString( ) ) )
                         {    to_delete.push_back( h.EdgeObjectIndexByIndexFrom( 
                                   v, j ) );    }    }    }
               if ( h.Sink(v) )
               {    for ( int j = 0; j < h.To(v).isize( ); j++ )
                    {    basevector b = kbb.Seq( h.EdgeObjectByIndexTo( v, j ) );
                         if ( !b.ToString( ).Contains( right.ToString( ) ) )
                         {    to_delete.push_back( h.EdgeObjectIndexByIndexTo( 
                                   v, j ) );    }    }    }    }
          UniqueSort(to_delete);
          h.DeleteEdges(to_delete);
          h.RemoveDeadEdgeObjects( );
          h.RemoveEdgelessVertices( );
          h.RemoveUnneededVertices( );
          if ( to_delete.empty( ) ) break;    }    }







// ---- computes a sorted vector of (kmer, kmer_pos) pairs

void build_kmer_pos_vec(const unsigned K,
			const BaseVec & bv, 
			vec<pair<unsigned, unsigned> > * kmer_pos_p)
{
  const unsigned nb = bv.size();
  if (nb >= K) {
    const unsigned nk = nb - K + 1;
    const unsigned kmer_mask = (1u << (K << 1)) - 1u;
    unsigned kmer = 0;
    
    // -- add k-1 first bases to kmer

    for (unsigned ib = 0; ib < K - 1; ib++)
      kmer = (kmer << 2) | bv[ib];

    // -- cycle through all read kmers and store them 

    for (unsigned ik = 0; ik < nk; ik++) {
      kmer = (kmer_mask & ((kmer << 2) | bv[ik + K - 1]));
      kmer_pos_p->push_back(make_pair(kmer, ik));
    }

    // -- sort for ease of future lookup

    sort(kmer_pos_p->begin(), kmer_pos_p->end());
  }
}







void PickPath( 
     // inputs:
     const int low, const int high, const efasta& etiglet, const basevector& tig, 
     const int estart, const int estop, const HyperKmerPath& h, 
     const KmerBaseBroker& kbb, const vecbasevector& basesy, 
     const vecqualvector& qualsy, const int maxread, const vec< vec<int> >& paths,
     // heuristics:
     const int max_expand_to, const int q_junk_max, const int max_non_losers,
     // logging:
     const int TIG, const Bool QLT, const String& tmp_dir, 
     const Bool SHOW_VOTES, const String& data_dir, ostream& rout,
     // output:
     vec< triple<int,int,String> >& replacements )
{
     int K = h.K( );
     vec<basevector> bpaths;
     rout << "\npaths:\n";
     for ( int m = 0; m < paths.isize( ); m++ )
     {    rout << "[" << m+1 << "] ";
          basevector b = kbb.Seq( h.EdgeObject( paths[m][0] ) );
          rout << "[";
          for ( int j = 0; j < paths[m].isize( ); j++ )
          {    if ( j > 0 ) rout << "-";
               rout << BaseAlpha( paths[m][j] );    
               if ( j > 0 )
               {    b.resize( b.isize( ) - (K-1) );
                    b = Cat( b, kbb.Seq( h.EdgeObject( paths[m][j] ) ) );    }    }
          bpaths.push_back(b);
          rout << "]\n";    }    
     
     // Now add in the assembly sequences.

     vec<basevector> A;
     if ( !etiglet.ExpandTo( A, max_expand_to ) ) rout << "expansion exploded\n";
     Bool first = True;
     for ( int j = 0; j < A.isize( ); j++ )
     {    if ( Member( bpaths, A[j] ) ) continue;
          if (first) rout << "extra paths from assembly:\n";
          first = False;
          rout << "[" << bpaths.size( ) + 1 << "] " << A[j] << "\n";
          bpaths.push_back( A[j] );    }
     if ( bpaths.empty( ) ) { rout << "no path\n"; return; }
     if ( bpaths.solo( ) ) { rout << "unique path\n"; return; }

     // Assess sequences versus reference.

     if (QLT)
     {    Mkdir777(tmp_dir);
          String qlt_fasta = tmp_dir + "/x.QLT." + "fasta";
          {    Ofstream( out, qlt_fasta );
               basevector x;
               x.Print(out); // dummy
               for ( int j = 0; j < bpaths.isize( ); j++ )
               {    const int big_flank = 8000;
                    int nleft = Min( big_flank, low );
                    int nright = Min( big_flank, tig.isize( ) - high );
                    basevector left( tig, low - nleft, nleft );
                    basevector right( tig, high, nright );
                    basevector b = Cat( left, bpaths[j], right );
                    b.Print( out, j+1 );    }    }
          SystemSucceed( "QueryLookupTable K=12 MM=12 MC=0.15 SEQS=" + qlt_fasta 
               + " L=" + data_dir + "/genome.lookup VISUAL=True PARSEABLE=True "
               "QUIET=True NH=True > " + tmp_dir + "/x.QLT.out" );
          vec<look_align> aligns;
          LoadLookAligns( tmp_dir + "/x.QLT.out", aligns );
          vec<int> errs( aligns.size( ) );
          for ( int j = 0; j < aligns.isize( ); j++ )
               errs[j] = aligns[j].Errors( );
          vec<int> ids( aligns.size( ), vec<int>::IDENTITY );
          SortSync( errs, ids );
          Bool winner_known = False;
          if ( errs.size( ) == 0 
               || ( errs.size( ) >= 2 && errs[0] == errs[1] )
               || !aligns[ ids[0] ].FullLength( ) )
          {    rout << "\nWarning: no clear winner\n";    }
          else
          {    winner_known = True;
               int w = ids[0];
               rout << "\nReference says winner is [" << w+1 
                    << "], errors = " << errs[0] << "\n";    }
          if ( !winner_known )
          {    rout << "\nAlignments:\n";
               fast_ifstream in( tmp_dir + "/x.QLT.out" );
               String line;
               while(1)
               {    getline( in, line );
                    if ( in.fail( ) ) break;
                    if ( !line.Contains( "QUERY", 0 ) ) 
                         rout << line << "\n";    }    }
          SystemSucceed( "/bin/rm -rf " + tmp_dir );    }

     // Extend paths using flanking sequence.

     vec<basevector> bpaths_orig(bpaths);
     int left_flank_size = Min( maxread, low );
     int right_flank_size = Min( maxread, tig.isize( ) - high );
     basevector left_flank( tig, low - left_flank_size, left_flank_size );
     basevector right_flank( tig, high, right_flank_size );
     for ( int j = 0; j < bpaths.isize( ); j++ )
          bpaths[j] = Cat( left_flank, bpaths[j], right_flank );

     // Align the reads to the extended paths, or "targets".  Note 
     // the following issues:
     // (1) We do not allow indels.
     // (2) We do not take into account consistency with the partner placement.
     //     Thus the best placement for the read taken individually might not be the
     //     best placement overall.
     // Voting.  For each gap-free placement of a read, we form the sum of its
     // quality scores at mismatches.  Quality scores of two are less are treated as
     // zero.  We find the target having the lowest score (the "winner"), and let
     // delta be the difference between the score of the runner up and the winner.
     // The target is then assigned a vote of 1 - 10^(-delta/10).

     rout << Date() << ": aligning" << endl;
     
     const size_t n_bpaths = bpaths.size();
     const size_t n_basesy = basesy.size();

     vec< vec<double> > VOTE(n_bpaths, vec<double>(n_bpaths, 0.0));
     vec< vec<unsigned> > SCORE(n_basesy, vec<unsigned>(n_bpaths, 1000000000));
     

     if (1) {
       const unsigned K = 7;
       const unsigned n_hits_min = 1;
       
       // ---- associate kmers in bpaths with their positions

       vec< vec<pair<unsigned, unsigned> > > kmer_pos_bpaths(n_bpaths);
       for (size_t i_bpaths = 0; i_bpaths < n_bpaths; i_bpaths++)
	 build_kmer_pos_vec(K, bpaths[i_bpaths], &(kmer_pos_bpaths[i_bpaths]));
       
       // ---- associate kmers in basesy with their positions

       vec< vec<pair<unsigned, unsigned> > > kmer_pos_basesy(n_basesy);
       for (size_t i_basesy = 0; i_basesy < n_basesy; i_basesy++)
	 build_kmer_pos_vec(K, basesy[i_basesy], &(kmer_pos_basesy[i_basesy]));
       
       // ---- cycle over all reads

       for (size_t i_basesy = 0; i_basesy < n_basesy; i_basesy++) {
	 const BaseVec & bv_r = basesy[i_basesy];
	 const QualVec & qv_r = qualsy[i_basesy];
	 const unsigned nb_r = bv_r.size();
	 const unsigned nk_r = nb_r - K + 1;
	 const vec<pair<unsigned, unsigned> > & kp_r = kmer_pos_basesy[i_basesy];

	 // ---- cycle over all paths
	 
	 for (size_t i_bpaths = 0; i_bpaths < n_bpaths; i_bpaths++) {
	   const BaseVec & bv_bp = bpaths[i_bpaths];
	   const unsigned nb_bp = bv_bp.size();
	   const unsigned nk_bp = nb_bp - K + 1;
	   const vec<pair<unsigned, unsigned> > & kp_bp = kmer_pos_bpaths[i_bpaths];
	   
	   // ---- associate with every path position the number of reads that 
	   //      align by a kmer.
	   //      a perfect read alignment at some path position would yield
	   //      a number of hits equal to the number of kmers in the read.

	   map<unsigned, unsigned> pos_hits;
	   for (unsigned ik_r = 0, ik0_bp = 0; ik_r < nk_r  &&  ik0_bp < nk_bp; ik_r++) {

	     const unsigned kmer_r = kp_r[ik_r].first;

	     unsigned ik_bp = ik0_bp;
	     while (ik_bp < nk_bp && kp_bp[ik_bp].first < kmer_r)
	       ik_bp++;

	     while (ik_bp < nk_bp && kp_bp[ik_bp].first == kmer_r) { // found some hits!
	       const unsigned & pos_bp = kp_bp[ik_bp].second;
	       const unsigned & pos_r  = kp_r[ik_r].second;
	       if (pos_bp >= pos_r && nk_bp - pos_bp >= nk_r - pos_r)
		 pos_hits[pos_bp - pos_r]++;
	       ik_bp++;
	     }

	     // next kmer is not the same: update ik0_bp

	     if (ik_r < nk_r - 1 && kp_r[ik_r + 1].first > kmer_r) 
	       ik0_bp = ik_bp;
	   }
       
	   // ---- compute the score for the read/path pair

	   unsigned & q_sum_min = SCORE[i_basesy][i_bpaths];


	   for (map<unsigned, unsigned>::iterator it = pos_hits.begin();
		it != pos_hits.end(); it++) {
	     if (it->second >= n_hits_min) {
	       unsigned ib0 = it->first;
	       unsigned q_sum = 0;
	       /*
	       cout << "i_bpaths= " << i_bpaths << " i_basesy= " << i_basesy << endl;
	       cout << "pos= " << it->first << " n_hits= " << it->second << " ib0= " << ib0 << endl;
	       cout << hieroglyphs(bv_r) << endl;
	       for (unsigned ib_r = 0; ib_r < nb_r; ib_r++)
                 cout << hieroglyph(bv_bp[ib0 + ib_r]);
	       cout << endl;
	       */

	       for (unsigned ib_r = 0; (ib_r < nb_r && q_sum < q_sum_min); ib_r++) {
                 if (bv_r[ib_r] != bv_bp[ib0 + ib_r] &&
		     qv_r[ib_r] > q_junk_max)
		   q_sum += qv_r[ib_r];
	       }
		      
	       if (q_sum < q_sum_min)
		 q_sum_min = q_sum;
	     }
	   }

		
	 }	 

       }


     }
     else {
       for (size_t i_basesy = 0; i_basesy < n_basesy; i_basesy++) {

	 const BaseVec & bv_r = basesy[i_basesy];
	 const QualVec & qv_r = qualsy[i_basesy];
	 const size_t nb_r = bv_r.size();
	 
	 for (size_t i_bpaths = 0; i_bpaths < n_bpaths; i_bpaths++) {
	 
	   const BaseVec & bv_bp = bpaths[i_bpaths];
	   const size_t nb_bp = bv_bp.size();

	   unsigned & q_sum_min = SCORE[i_basesy][i_bpaths];

	   for (unsigned ib = 0; ib <= nb_bp - nb_r; ib++) {
	     unsigned q_sum = 0;
	   
	     for (unsigned ib_r = 0; (ib_r < nb_r && q_sum < q_sum_min); ib_r++) {
               if (bv_r[ib_r] != bv_bp[ib + ib_r] &&
	           qv_r[ib_r] > q_junk_max)
                 q_sum += qv_r[ib_r];
	     }
	   
	     if (q_sum < q_sum_min)
	       q_sum_min = q_sum;
		
	   }
         }
       }

     }

     rout << Date() << ": voting" << endl;
       

     for (unsigned i1_bpaths = 0; i1_bpaths < n_bpaths; i1_bpaths++) {
       for (unsigned i2_bpaths = i1_bpaths + 1; i2_bpaths < n_bpaths; i2_bpaths++) {
	 
	 double & vote12 = VOTE[i1_bpaths][i2_bpaths];
	 double & vote21 = VOTE[i2_bpaths][i1_bpaths];
	 vote12 = vote21 = 0;

	 for (unsigned i_basesy = 0; i_basesy < n_basesy; i_basesy++) {

	   const unsigned & score1 = SCORE[i_basesy][i1_bpaths];
	   const unsigned & score2 = SCORE[i_basesy][i2_bpaths];
	   
	   if      (score1 > score2) vote21 += 1.0 - pow(10.0, -0.1 * (score1 - score2));
	   else if (score2 > score1) vote12 += 1.0 - pow(10.0, -0.1 * (score2 - score1));
	 }
       }
     }


     rout << Date() << ": finding winner" << endl;
     






     
     int winner = -1;
     for ( int i1 = 0; i1 < bpaths.isize( ); i1++ )
     {    Bool i1_wins = True;
          for ( int i2 = 0; i2 < bpaths.isize( ); i2++ )
          {    if ( i1 == i2 ) continue;
               if ( VOTE[i1][i2] == 0 || VOTE[i2][i1] >= 2.0 
                    || VOTE[i1][i2] < 10.0 * VOTE[i2][i1] )
               {    i1_wins = False;    }    }
          if (i1_wins) winner = i1;    }
     if ( winner >= 0 )
     {    rout << "winner: [" << winner+1 << "]\n";
          if ( bpaths_orig[winner].ToString( ) == etiglet )
               rout << "equals assembly\n";
          else
          {    rout << "replacing assembly sequence\n" << etiglet 
                    << "\nby\n" << bpaths_orig[winner].ToString( ) << "\n";
               PRINT3_TO( rout, TIG, estart, estop );
               #pragma omp critical
               {    replacements.push( estart, estop,
                         bpaths_orig[winner].ToString( ) );    }    }
          return;    }

     vec<Bool> loser( bpaths.size( ), False );
     vec<Bool> beats_something( bpaths.size( ), False );
     vec<Bool> something_beats( bpaths.size( ), False );

     for ( int i1 = 0; i1 < bpaths.isize( ); i1++ )
     {    for ( int i2 = 0; i2 < bpaths.isize( ); i2++ )
          {    if ( i1 == i2 ) continue;
               if ( VOTE[i2][i1] < 2.0 && VOTE[i1][i2] >= 10.0 * VOTE[i2][i1] )
               {    beats_something[i1] = True;
                    something_beats[i2] = True;    }    }    }
     for ( int i1 = 0; i1 < bpaths.isize( ); i1++ )
          if ( something_beats[i1] && !beats_something[i1] ) loser[i1] = True;
     rout << "non-losers:\n";
     int n_non_losers = 0;
     for ( int i1 = 0; i1 < bpaths.isize( ); i1++ )
     {    if ( !loser[i1] ) 
          {    rout << "[" << i1+1 << "]\n";
               n_non_losers++;    }    }
     for ( int i1 = 0; i1 < bpaths.isize( ); i1++ )
     for ( int i2 = i1 + 1; i2 < bpaths.isize( ); i2++ )
     {    if ( !loser[i1] && !loser[i2] )
          {    rout << "vote:\n" << "[" << i1+1 << "] " << VOTE[i1][i2] << "\n";
               rout << "[" << i2+1 << "] " << VOTE[i2][i1] << "\n";    }    }
     if ( n_non_losers <= max_non_losers )
     {    vec<basevector> nl;
          for ( int j = 0; j < bpaths.isize( ); j++ )
               if ( !loser[j] ) nl.push_back( bpaths_orig[j] );
          Sort(nl), Sort(A);
          if ( nl != A )
          {    rout << "replacing assembly sequence\n" << etiglet 
                    << "\nby\n" << efasta(nl) << "\n";
               PRINT3_TO( rout, TIG, estart, estop );
               #pragma omp critical
               {    replacements.push( estart, estop,
                         efasta(nl) );    }    }    }    }

class heuristics {

     public:

     int min_calls;
     double agree_ceil_weak;
     double agree_floor_strong;
     double MIN_FRAC;
     int q_junk_max;
     int max_tandem_period;
     int min_tandem_copies;
     int min_tandem_length;
     double max_offby;
     int max_non_losers;
     int max_expand_to;
     int max_paths;
     double gapdev_mult;
     int boundary_push;
     int max_cyclic_paths;
     int max_cyclic_loops;
     double min_devs_off;
     Bool do_cyclic;
     int max_cyclic_edges;

};

class logging_control {

     public:

     String DUMP_EC;
     Bool SHOW_VOTES;
     Bool QLT;
     String data_dir;
     String DOT;

};

void ProcessWindow( 
     // inputs:
     const heuristics& heur, const logging_control& logc, const String& tmp_dir,
     pair<int,int>& window, const String& signal,
     const vec< triple<ho_interval,int,int> >& fragbounds, const basevector& tig,
     const efasta& T, const int TIG, const int K, const vec<int>& ids_frag,
     const vec<int>& ids_jump, const vec<read_loc>& locs,
     const vec<int>& loc_id_frag, const vec<int>& loc_id_jump,
     const vecbasevector& bases, const vecqualvector& quals,
     const vec<int>& readlengths,
     // outputs:
     ostream& rout, vec< triple<int,int,int> >& markups,
     vec< triple<int,int,String> >& replacements )
{
     rout << signal << "\n";
     restart:
     int low = window.first, high = window.second, nbounds = 0;
     vec<int> gap, gapdev;
     for ( int j = 0; j < fragbounds.isize( ); j++ )
     {    if ( Subset( ho_interval( low, high ), fragbounds[j].first ) ) 
          {    nbounds++;
               gap.push_back( fragbounds[j].second );
               gapdev.push_back( fragbounds[j].third );    }    }
     int gap_ave, gapdev_ave;
     GapStatsAlt( gap, gapdev, gap_ave, gapdev_ave );
     PRINT3_TO( rout, nbounds, gap_ave, gapdev_ave );
     double devs_off = -1;
     if ( nbounds > 0 ) devs_off = Abs( double(gap_ave)/double(gapdev_ave) );
     Bool way_off = ( devs_off >= heur.min_devs_off );
     rout << "\n";
 
     // Get the contig chunk and its flanks.

     basevector left( tig, low, K ), right( tig, high - K, K );
     basevector tiglet( tig, low, high - low );
     int pos = 0, estart = -1, estop = -1;
     for ( int j = 0; j < T.isize( ); j++ )
     {    if ( pos == low ) estart = j;
          if ( pos == high ) estop = j;
          if ( T[j] == '{' )
          {    for ( j++; j < T.isize( ); j++ )
               {    if ( T[j] == ',' ) break;
                    pos++;    }
               for ( j++; j < T.isize( ); j++ )
                    if ( T[j] == '}' ) break;    }
          else pos++;    }
     efasta etiglet = efasta( T.substr( estart, estop - estart ) );
     etiglet.Print( rout, "etiglet (" + ToString( etiglet.Length1( ) ) + " bases)" );

     // Isolate the reads for this region.

     vecbasevector basesy;
     vecqualvector qualsy;
     for ( int j = 0; j < ids_frag.isize( ); j++ )
     {    const read_loc& rl = locs[ loc_id_frag[j] ];
          if ( IntervalOverlap( rl.Start( ), rl.Stop( ), low, high ) > 0 )
          {    basesy.push_back( bases[j] );
               qualsy.push_back( quals[j] );    }    }
     for ( int j = 0; j < ids_jump.isize( ); j++ )
     {    const read_loc& rl = locs[ loc_id_jump[j] ];
          if ( IntervalOverlap( rl.Start( ), rl.Stop( ), low, high ) > 0 )
          {    basesy.push_back( bases[ j + ids_frag.isize( ) ] );
               qualsy.push_back( 
                    quals[ j + ids_frag.isize( ) ] );    }    }

     // Correct errors in the reads.

     BaseVecVec basesx = basesy;
     {    const size_t nqv = qualsy.size();
          QualNibbleVecVec qualsx(nqv);
          for (size_t iqv = 0; iqv < nqv; iqv++)
          {    for (size_t iq = 0; iq < qualsy[iqv].size(); iq++)
               {    qualsx[iqv].push_back(qualsy[iqv][iq]);    }    }
          const unsigned K = 24;
          const unsigned n_cycles = 2;
          const unsigned verbosity = 0;
          find_errors(efp_default, K, &basesx, &qualsx, n_cycles, verbosity);
          if ( logc.DUMP_EC != "" ) basesx.WriteAll(logc.DUMP_EC);    }

     vecKmerPath Paths, Paths_rc, unipaths;
     vec<tagged_rpint> Pathsdb, unipathsdb;
     ReadsToPathsCoreY( basesx, K, Paths );
     CreateDatabase( Paths, Paths_rc, Pathsdb );
     Unipath( Paths, Paths_rc, Pathsdb, unipaths, unipathsdb );
     KmerBaseBroker kbb( K, Paths, Paths_rc, Pathsdb, basesx );
     digraph AA;
     BuildUnipathAdjacencyGraph(Paths, Paths_rc, Pathsdb, unipaths, unipathsdb, AA);
     HyperKmerPath h;
     BuildUnipathAdjacencyHyperKmerPath( K, AA, unipaths, h );
     int MIN_COMPONENT = 20;
     if ( MIN_COMPONENT > 0 ) h.RemoveSmallComponents(MIN_COMPONENT);
     h.RemoveDeadEdgeObjects( );
     h.RemoveEdgelessVertices( );
     h.RemoveUnneededVertices( );

     // Delete initial edges that do not contain the left flank, and
     // terminal edges that do no contain the right flank.

     DeleteIllegalTerminators( h, kbb, left, right );

     // If an initial [resp. terminal] edge contains the left [resp. right]
     // flank, truncate it there.

     for ( int v = 0; v < h.N( ); v++ )
     {    if ( h.Source(v) )
          {    for ( int j = 0; j < h.From(v).isize( ); j++ )
               {    KmerPath& e = h.EdgeObjectByIndexFromMutable( v, j );
                    basevector b = kbb.Seq(e);
                    int p = b.ToString( ).Position( left.ToString( ) );
                    if ( p > 0 )
                    {    KmerPath enew;
                         e.CopySubpath( e.Begin( ) + p, e.End( ), enew );
                         e = enew;    }    }    }
          if ( h.Sink(v) )
          {    for ( int j = 0; j < h.To(v).isize( ); j++ )
               {    KmerPath& e = h.EdgeObjectByIndexToMutable( v, j );
                    basevector b = kbb.Seq(e);
                    int p = b.ToString( ).Position( right.ToString( ) );
                    if ( p >= 0 && p + K < b.isize( ) )
                    {    KmerPath enew;
                         e.CopySubpath( e.Begin( ), 
                              e.End( ) - ( b.isize( ) - p - K), enew );
                         e = enew;    }    }    }    }

     // Prune branches using MIN_FRAC.

     vec<double> cov;
     Coverage( h, Pathsdb, cov );
     vec<int> to_delete;
     for ( int v = 0; v < h.N( ); v++ )
     {    for ( int j1 = 0; j1 < h.From(v).isize( ); j1++ )
          {    for ( int j2 = 0; j2 < h.From(v).isize( ); j2++ )
               {    int e1 = h.EdgeObjectIndexByIndexFrom( v, j1 );
                    int e2 = h.EdgeObjectIndexByIndexFrom( v, j2 );
                    // if ( h.From(v)[j1] != h.From(v)[j2] ) continue;
                    double c1 = cov[e1], c2 = cov[e2];
                    if ( c1 >= 2.0 && c2 >= 2.0 ) continue;
                    if ( c1 < int( floor( heur.MIN_FRAC * (c1+c2) ) ) )
                         to_delete.push_back(e1);    }    }
          for ( int j1 = 0; j1 < h.To(v).isize( ); j1++ )
          {    for ( int j2 = 0; j2 < h.To(v).isize( ); j2++ )
               {    int e1 = h.EdgeObjectIndexByIndexTo( v, j1 );
                    int e2 = h.EdgeObjectIndexByIndexTo( v, j2 );
                    // if ( h.To(v)[j1] != h.To(v)[j2] ) continue;
                    double c1 = cov[e1], c2 = cov[e2];
                    if ( c1 >= 2.0 && c2 >= 2.0 ) continue;
                    if ( c1 < int( floor( heur.MIN_FRAC * (c1+c2) ) ) )
                         to_delete.push_back(e1);    }    }    }
     UniqueSort(to_delete);
     h.DeleteEdges(to_delete);
     if ( MIN_COMPONENT > 0 ) h.RemoveSmallComponents(MIN_COMPONENT);
     h.RemoveDeadEdgeObjects( );
     h.RemoveEdgelessVertices( );
     h.RemoveUnneededVertices( );
     h.RemoveDeadEdgeObjects( );

     // Look for grubby kmers.

     vec<int> grubby;
     vec<basevector> B, BX;
     vec<qualvector> Q;
     for ( size_t j = 0; j < basesy.size( ); j++ )
     {    for ( int l = 0; l <= basesy[j].isize( ) - K; l++ )
          {    basevector b( basesy[j], l, K );
               qualvector q;
               q.SetToSubOf( qualsy[j], l, K );
               B.push_back(b), Q.push_back(q);    }    }
     for ( int j = 0; j < h.EdgeObjectCount( ); j++ )
     {    basevector b = kbb.Seq( h.EdgeObject(j) );
          for ( int l = 0; l <= b.isize( ) - K; l++ )
               BX.push( b, l, K );    }
     SortSync( B, Q ), Sort(BX);
     for ( int j = 0; j < h.EdgeObjectCount( ); j++ )
     {    basevector b = kbb.Seq( h.EdgeObject(j) );
          for ( int l = 0; l <= b.isize( ) - K; l++ )
          {    basevector c( b, l, K );
               for ( int r = 0; r < K; r++ )
               {    int m0 = c[r];
                    vec<Bool> ref( 4, False );
                    for ( int m = 0; m < 4; m++ )
                    {    if ( m == m0 ) ref[m] = True;
                         else 
                         {    basevector d(c);
                              d.Set( r, m );
                              if ( BinMember( BX, d ) ) ref[m] = True;    }    }
                    if ( Sum(ref) == 1 ) continue;
                    vec<int> low(4), high(4);
                    for ( int m = 0; m < 4; m++ )
                    {    basevector d(c);
                         d.Set( r, m );
                         low[m] = lower_bound( B.begin( ), B.end( ), d )
                              - B.begin( );
                         high[m] = upper_bound( B.begin( ), B.end( ), d )
                              - B.begin( );    }
                    Bool have_nonref = False;
                    for ( int m = 0; m < 4; m++ )
                    {    if ( m == m0 ) continue;
                         if ( low[m] < high[m] ) have_nonref = True;    }
                    if ( !have_nonref ) continue;

                    // Kill the kmer if it is a mutation of another kmer,
                    // if the other kmer has at least 10 times as much 
                    // support, and if the maximum quality score at the 
                    // mutation position is < the median quality score of 
                    // the other kmer.
     
                    const int cov_qtest_ratio = 10;
                    vec< vec<int> > qs(4);
                    for ( int m = 0; m < 4; m++ )
                    {    for ( int x = low[m]; x < high[m]; x++ )
                              qs[m].push_back( Q[x][r] );
                         Sort( qs[m] );    }
                    for ( int m = 0; m < 4; m++ )
                    {    if ( qs[m0].size( ) * cov_qtest_ratio 
                              > qs[m].size( ) )
                         {    continue;    }
                         if ( qs[m0].nonempty( ) && qs[m0].back( ) <
                              Median( qs[m] ) )
                         {    grubby.push_back(j);    }    }
                    Bool ec_verbose = False;
                    if (ec_verbose)
                    {    rout << "\n" << c << " " << r << "\n";
                         for ( int m = 0; m < 4; m++ )
                         {    const vec<int>& v = qs[m];
                              rout << as_base(m) << (m == c[r] ? " (ref" : "");
                              if ( v.isize( ) < 20 )
                              {    for ( int x = 0; x < v.isize( ); x++ )
                                        rout << " " << v[x];    }
                              else
                              {    rout << " " << v.size( ) << " quals of median " 
                                        << v[ v.size( )/2 ];    }
                              rout << "\n";    }    }    }    }    }
     UniqueSort(grubby);
     h.DeleteEdges(grubby);
     if ( MIN_COMPONENT > 0 ) h.RemoveSmallComponents(MIN_COMPONENT);
     h.RemoveDeadEdgeObjects( );
     h.RemoveEdgelessVertices( );
     h.RemoveUnneededVertices( );
     h.RemoveDeadEdgeObjects( );

     // Delete initial edges that do not contain the left flank, and
     // terminal edges that do no contain the right flank.

     DeleteIllegalTerminators( h, kbb, left, right );

     // If left and right are not both present, kill the graph.

     int left_id = -1, right_id = -1;
     vec<int> to_left, to_right;
     h.ToLeft(to_left), h.ToRight(to_right);
     for ( int j = 0; j < h.EdgeObjectCount( ); j++ )
     {    if ( kbb.Seq( h.EdgeObject(j) ).ToString( ).Contains( left.ToString( ) ) )
          {    left_id = to_left[j];    }
          if ( kbb.Seq( h.EdgeObject(j) ).ToString( ).Contains( right.ToString( ) ) )
          {    right_id = to_right[j];    }    }
     if ( left_id < 0 || right_id < 0 )
     {    rout << "graph is empty" << endl;
          if ( !etiglet.Contains( "{" ) )
          {    rout << "marking assembly sequence\n";
               PRINT3_TO( rout, TIG, estart, estop );
               #pragma omp critical
               {    markups.push( TIG, estart, estop );    }    }
          return;    }

     // Present results.

     rout << "\n";
     PRINT_TO( rout, h.EdgeObjectCount( ) );
     if ( logc.DOT != "" )
     {    Ofstream( out, logc.DOT );
          h.PrintSummaryDOT0w( out, False, False, True, NULL, False );    }
     Coverage( h, Pathsdb, cov );
     for ( int e = 0; e < h.EdgeObjectCount( ); e++ )
     {    int v = to_left[e], w = to_right[e];
          rout << ">edge_" << BaseAlpha(e) << " " << v << ":" << w 
               << " kmers=" << h.EdgeLengthKmers(e)
               << " cov=" << ToString( cov[e] ) << "\n";
          kbb.ToSequence( h.EdgeObject(e) ).PrintN(rout);    }    

     // Handle the case where there are cycles.

     if ( !h.Acyclic( ) )
     {    rout << "cyclic\n";

          // Give up if too many edges.

          if ( h.EdgeObjectCount( ) > heur.max_cyclic_edges )
          {    rout << "\ntoo many edges in cyclic graph, giving up\n";
               return;    }

          // To proceed we require that the graph has a unique source and
          // sink, and unique edges emanating from them.  If we don't
          // find a source or a sink, we extend the boundaries of the
          // window, which is done via a goto.  Yecch.

          vec<int> sources, sinks;
          h.Sources(sources), h.Sinks(sinks);
          PRINT2_TO( rout, sources.size( ), sinks.size( ) );
          if ( sources.empty( ) )
          {    if ( low == 0 )
               {    rout << "\nno sources, can't push left boundary "
                         << "further, giving up\n";
                    return;    }
               low -= heur.boundary_push;
               if ( low < 0 ) low = 0;
               window.first = low;
               rout << "\n*** no sources, pushing left boundary of "
                    << "window to " << low << " ***\n\n";
               goto restart;    }
          if ( sinks.empty( ) )
          {    if ( high == tig.isize( ) )
               {    rout << "\nno sinks, can't push right boundary "
                         << "further, giving up\n";
                    return;    }
               high += heur.boundary_push;
               if ( high > tig.isize( ) ) high = tig.size( );
               window.second = high;
               rout << "\n*** no sinks, pushing right boundary of window "
                    << "to " << high << " ***\n\n";
               goto restart;    }
          if ( !sources.solo( ) || !sinks.solo( ) )
          {    rout << "don't have unique source and sink\n";
               return;    }
          int v = sources[0], w = sinks[0];
          int source_edges = h.From(v).size( );
          int sink_edges = h.To(w).size( );
          if ( source_edges != 1 || sink_edges != 1 )
          {    rout << "have unique source and sink\n";
               PRINT2_TO( rout, source_edges, sink_edges );
               rout << "don't have unique edges emanating from "
                    << "source and sink\n";
               return;    }
          rout << "have unique source and sink and edges "
               << "emanating from them\n";
          rout << BaseAlpha( h.EdgeObjectIndexByIndexFrom( v, 0 ) )
               << " --> ... --> "
               << BaseAlpha( h.EdgeObjectIndexByIndexTo( w, 0 ) ) << "\n";
          HyperBasevector hb( h, kbb );
          int in_id = hb.EdgeObjectIndexByIndexFrom( v, 0 );
          int out_id = hb.EdgeObjectIndexByIndexTo( w, 0 );
          if ( !hb.EdgeObject(in_id).ToString( ).Contains( 
               left.ToString( ), 0 )
               || !hb.EdgeObject(out_id).ToString( ).Contains( 
               right.ToString( ), -1 ) )
          {    rout << "can't find left and right\n";
               return;    }

          // Could there be only a small number of paths through the
          // graph that are consistent with distance constraints?

          if ( nbounds > 0 )
          {    vec<int> L;
               for ( int j = 0; j < h.EdgeObjectCount( ); j++ )
                    L.push_back( h.EdgeLengthKmers(j) );
               digraphE<int> G( h, L );
               int ex = etiglet.Length1( ) - (K-1);
               int L1 = Max( 0, ex - int( ceil( double(gap_ave) 
                    + heur.gapdev_mult * double(gapdev_ave) ) ) );
               int L2 = Max( 0, ex - int( floor( double(gap_ave) 
                    - heur.gapdev_mult * double(gapdev_ave) ) ) );
               vec< vec<int> > paths;
               if ( G.AllPathsLengthRange( v, w, L1, L2, to_right, paths,
                    heur.max_cyclic_paths, heur.max_cyclic_loops ) )
               {    rout << "\nsee " << paths.size( ) 
                         << " paths through the graph that are "
                         << "consistent with distance constraints\n";    
                    if ( paths.nonempty( ) )
                    {    PickPath( low, high, etiglet, tig, estart, estop, h, kbb, 
                              basesy, qualsy, Max(readlengths), paths, 
                              heur.max_expand_to, heur.q_junk_max, 
                              heur.max_non_losers, 
                              TIG, logc.QLT, tmp_dir, logc.SHOW_VOTES, 
                              logc.data_dir, rout, replacements );    
                         return;    }    }    }

          // Test for simple cyclic graph.

          Bool simple_cyclic = False;
          if ( h.EdgeObjectCount( ) == 4 )
          {    if ( h.From(v).solo( ) && h.To(w).solo( ) )
               {    int x = h.From(v)[0], y = h.To(w)[0];
                    if ( h.From(x).solo( ) && h.To(x).size( ) == 2
                         && h.From(y).size( ) ==2 && h.To(y).solo( ) )
                    {    simple_cyclic = True;
                         rout << "simple cyclic\n";    }    }    }
          if ( !simple_cyclic && !heur.do_cyclic ) 
          {    if (way_off) 
               {    rout << "way off\n" << "marking assembly sequence\n";
                    PRINT3_TO( rout, TIG, estart, estop );
                    #pragma omp critical
                    {    markups.push( TIG, estart, estop );    }    }
               rout << "giving up, as the graph is cyclic but not simple" << endl;
               return;    }
     
          // Create an extended version of the graph in which we add
          // a readlength of bases to both ends.

          HyperBasevector hbplus(hb);
          basevector& in_edge = hbplus.EdgeObjectMutable(in_id);
          basevector& out_edge = hbplus.EdgeObjectMutable(out_id);
          int mr = Max(readlengths);
          int lext = Min( mr, low ), rext = Min( mr, tig.isize( ) - high );
          basevector left_ext( tig, low - lext, lext );
          basevector right_ext( tig, high, rext );
          in_edge = Cat( left_ext, in_edge );
          out_edge = Cat( out_edge, right_ext );

          // Align the reads to the edges.

          int source_sinks = 0;
          Bool exploded = False;
          vec< vec< vec<int> > > subpaths( readlengths.size( ) );
          for ( int r = 0; r < readlengths.isize( ); r++ )
          {    exploded = !AllSubpathTargets( 
                    hbplus, readlengths[r], subpaths[r] );
               if (exploded)
                    rout << "Warning, AllSubpathTargets exploded\n";    }
          for ( size_t r = 0; r < basesy.size( ); r++ )
          {    const basevector& R = basesy[r];
               const qualvector& Q = qualsy[r];
               int px = BinPosition( readlengths, R.isize( ) );
               vec<int> score( subpaths[px].size( ) );
               for ( int j = 0; j < subpaths[px].isize( ); j++ )
               {    const vec<int>& p = subpaths[px][j];
                    basevector target = hbplus.EdgeObject( p[0] );
                    for ( int l = 1; l < p.isize( ); l++ )
                    {    target.resize( target.isize( ) - (K-1) );
                         target = Cat( 
                              target, hbplus.EdgeObject( p[l] ) );    }
                    int n1 = hbplus.EdgeLengthBases( p[0] );
                    int bestq = 1000000000;
                    int nx = n1;
                    for ( int l = 1; l < p.isize( ) - 1; l++ )
                         nx += hbplus.EdgeLengthKmers( p[l] );
                    for ( int start = Max( 0, nx - R.isize( ) + 1 ); 
                         start < n1 - (K-1); start++ )
                    {    if ( start + R.isize( ) - 1 < target.isize( ) )
                         {    int qsum = 0;
                              for ( int l = 0; l < R.isize( ); l++ )
                              {    if ( R[l] != target[start+l] ) 
                                        qsum += Q[l];    }
                              bestq = Min( bestq, qsum );    }    }
                    score[j] = bestq;    }    
               vec<int> ids( subpaths[px].size( ), vec<int>::IDENTITY );
               SortSync( score, ids );
               const int max_diff = 20;
               const int max_score = 200;
               if ( score[0] <= max_score )
               {    rout << "\nalignments of read " << r << "\n";
                    rout << R.ToString( ) << "\n";
                    Bool source_sink = True;
                    vec< vec<String> > rows;
                    for ( int j = 0; j < score.isize( ); j++ )
                    {    if ( score[j] > score[0] + max_diff ) break;
                         vec<String> row;
                         row.push_back( "[" + ToString(ids[j]+1) + "]" );
                         row.push_back( 
                              "(score=" + ToString(score[j]) + ")" );   
                         const vec<int>& p = subpaths[px][ ids[j] ];
                         String r;
                         for ( int l = 0; l < p.isize( ); l++ )
                         {    r += ( l > 0 ? "-" : "" ) 
                                   + BaseAlpha( p[l] );    }
                         row.push_back(r);
                         rows.push_back(row);
                         if ( !h.Source( to_left[ p.front( ) ] )
                              || !h.Sink( to_right[ p.back( ) ] ) )
                         {    source_sink = False;    }    }
                    PrintTabular( rout, rows, 1, "lll" );
                    if (source_sink)
                    {    source_sinks++;
                         rout << "This read aligns as a "
                              << "source-sink.\n";    }    }    }
          rout << "\n";
          PRINT_TO( rout, source_sinks );
          const int min_source_sinks = 3;
          if ( source_sinks < min_source_sinks )
          {    rout << "Can't find enough source-sink edges, suggest "
                    << "using gap_ave = " << gap_ave << " +/- " 
                    << gapdev_ave 
                    << " to tweak\nnumber of copies of repeat.\n";
               if ( simple_cyclic && !exploded && nbounds > 0
                    && !signal.Contains( "ambiguity" ) )
               {    rout << "to perturb tandem count\n";    
                              
                    // Unwind etiglet along graph.

                    int x = h.From(v)[0];
                    int fwd_id = h.EdgeObjectIndexByIndexFrom( x, 0 );
                    int y = h.From(x)[0];
                    int rev_id = -1;
                    for ( int j = 0; j < h.From(y).isize( ); j++ )
                    {    if ( h.From(y)[j] == x )
                         {    rev_id = h.EdgeObjectIndexByIndexFrom( 
                                   y, j );    }    }
                    int xcount = -1;
                    for ( int count = 0; ; count++ )
                    {    vec<int> v;
                         v.push_back(in_id);
                         for ( int r = 0; r < count; r++ )
                              v.push_back( fwd_id, rev_id );
                         v.push_back( fwd_id, out_id );
                         basevector e = hb.EdgePathToBases(v);
                         if ( e.isize( ) > etiglet.isize( ) ) break;
                         if ( e.ToString( ) == etiglet )
                         {    xcount = count;
                              break;    }    }
                    if ( xcount < 0 )
                    {    rout << "unable to convert etiglet\n";
                         return;   }
                    int repeat_size = h.EdgeLengthKmers(fwd_id)
                         + h.EdgeLengthKmers(rev_id);
                    rout << "observe " << xcount << " repeat copies, "
                         << "each of size " << repeat_size << "\n";
                    int add_low = int( floor( ( double(-gap_ave) 
                         - heur.gapdev_mult * double(gapdev_ave) )
                         / double(repeat_size) ) );
                    int add_high = int( ceil( ( double(-gap_ave) 
                         + heur.gapdev_mult * double(gapdev_ave) )
                         / double(repeat_size) ) );
                    rout << "to add between " << add_low << " and "
                         << add_high << " copies of repeat\n";    
                    int new_count_low = Max( 0, xcount + add_low );
                    int new_count_high = Max( 0, xcount + add_high );
                    if ( !( new_count_low <= new_count_high ) )
                    {    rout << "doesn't make sense\n";
                         return;    }

                    // Replace tandem repeat.

                    vec<basevector> answer;
                    for ( int c = new_count_low; c <= new_count_high; c++ )
                    {    vec<int> v;
                         v.push_back(in_id);
                         for ( int r = 0; r < c; r++ )
                              v.push_back( fwd_id, rev_id );
                         v.push_back( fwd_id, out_id );
                         answer.push_back( hb.EdgePathToBases(v) );    }
                    rout << "replacing assembly sequence\n" 
                         << etiglet << "\nby\n" << efasta(answer) << "\n";
                    PRINT3_TO( rout, TIG, estart, estop );
                    #pragma omp critical
                    {    replacements.push( estart, estop,
                              efasta(answer) );    }    }    }
          return;    }

     // Handle the acyclic case.  First find all paths.

     rout << "acyclic" << endl;
     vec< vec<int> > paths;
     if ( !h.EdgePaths( left_id, right_id, paths, -1, heur.max_paths ) )
          rout << "path count too large\n";
     else
     {    PickPath( low, high, etiglet, tig, estart, estop, h, kbb, basesy, qualsy, 
               Max(readlengths), paths, heur.max_expand_to, heur.q_junk_max, 
               heur.max_non_losers, TIG, logc.QLT, tmp_dir, logc.SHOW_VOTES, 
               logc.data_dir, rout, replacements );    }    }

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(SUBDIR, "test");
     CommandArgument_String(SCAFFOLDS_IN);
     CommandArgument_String_OrDefault(SCAFFOLDS_OUT, SCAFFOLDS_IN + ".local");
     CommandArgument_Int_OrDefault_Doc(NUM_THREADS, -1,
          "number of threads to be used (use all available if negative)");
     CommandArgument_Int_OrDefault_Doc(VERBOSITY, 0, "can be 0, 1, or 2");
     CommandArgument_Bool_OrDefault(QLT, False);
     CommandArgument_Bool_OrDefault(WRITE, True);
     CommandArgument_Int_OrDefault_Doc(TIGX, -1,
          "if specified, process just this contig");
     CommandArgument_Int_OrDefault(STARTX, -1);
     CommandArgument_Int_OrDefault(STOPX, -1);
     CommandArgument_Bool_OrDefault_Doc(DIRECT, False,
          "log directly to cout, forces VERBOSITY >= 1 and NUM_THREADS = 1");
     CommandArgument_Bool_OrDefault(SHOW_VOTES, False);
     CommandArgument_Int_OrDefault(MAX_ASSEMBLY, 10000000);
     CommandArgument_Bool_OrDefault_Doc(RUN_ON_DIPLOID, False,
          "run on diploid genomes; default is to exit on diploids");
     CommandArgument_Bool_OrDefault(DO_CYCLIC, False);
     CommandArgument_Int_OrDefault(WINDOW_START, -1);
     CommandArgument_String_OrDefault_Doc(DUMP_EC, "",
          "file to dump error-corrected reads to, only recommended if used with "
          "TIGX+WINDOW_START so there is only one window analyzed");
     CommandArgument_String_OrDefault_Doc(DOT, "",
          "if specified, name of file to dump DOT graph to; should only be\n"
          "run with a single window");
     EndCommandArguments;

     // Define directories, etc.

     double clock = WallClockTime( );
     String data_dir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
     String head = sub_dir + "/" + SCAFFOLDS_IN;
     if ( STARTX >= 0 || STOPX >= 0 ) ForceAssertLt( STARTX, STOPX );
     if ( WINDOW_START >= 0 )
     {    ForceAssertGe( TIGX, 0 );
          const int default_wing = 6000;
          if ( STARTX < 0 ) STARTX = WINDOW_START - default_wing;
          if ( STOPX < 0 ) STOPX = WINDOW_START + default_wing;    }
     if ( TIGX >= 0 ) DIRECT = True;

     // Thread control

     if (DIRECT) NUM_THREADS = 1;
     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );

     // Define heuristics.

     const int K = 20;
     heuristics heur;
     heur.min_calls = 10;
     heur.agree_ceil_weak = 0.5;
     heur.agree_floor_strong = 0.9;
     heur.MIN_FRAC = 0.1;
     heur.q_junk_max = 2;
     heur.max_tandem_period = 8;
     heur.min_tandem_copies = 3;
     heur.min_tandem_length = 12;
     heur.max_offby = 10.0;
     heur.max_non_losers = 3;
     heur.max_expand_to = 2000;
     heur.max_paths = 500;
     heur.gapdev_mult = 3.0;
     heur.boundary_push = 20;
     heur.max_cyclic_paths = 5;
     heur.max_cyclic_loops = 100;
     heur.min_devs_off = 10.0;
     heur.do_cyclic = DO_CYCLIC;
     heur.max_cyclic_edges = 50;

     // Define logging control.

     logging_control logc;
     logc.SHOW_VOTES = SHOW_VOTES;
     logc.DUMP_EC = DUMP_EC;
     logc.QLT = QLT;
     logc.data_dir = data_dir;
     logc.DOT = DOT;

     // Load assembly.

     String efasta_file = sub_dir + "/" + SCAFFOLDS_IN + ".contigs.efasta";
     VecEFasta tigse;
     LoadEfastaIntoStrings( efasta_file, tigse );
     vecbasevector tigs( sub_dir + "/" + SCAFFOLDS_IN + ".contigs.fastb" );
     String supers_file = sub_dir + "/" + SCAFFOLDS_IN + ".superb";
     vec<superb> scaffolds;
     ReadSuperbs( supers_file, scaffolds );
     const int ploidy = FirstLineOfFile( data_dir + "/ploidy" ).Int( );

     // If assembly is too large, quit.

     int64_t total_size = 0;
     for ( size_t i = 0; i < tigs.size( ); i++ )
          total_size += tigs[i].size( );
     if ( total_size > MAX_ASSEMBLY || ( !RUN_ON_DIPLOID && ploidy > 1 ) )
     {    cout << "Assembly size exceeds MAX_ASSEMBLY, or ploidy > 1, "
               << "making no changes." << endl;
          Assembly A( scaffolds, tigse );
          A.WriteAll( sub_dir + "/" + SCAFFOLDS_OUT );
          SystemSucceed( "touch " + sub_dir + "/" + SCAFFOLDS_OUT + ".markup" );
          return 0;    }

     // Set up to access read locations.

     read_locs_on_disk locs_file( head, run_dir );

     // Compute distributions (turned off for now).

     vec<ProbFuncIntDist> pfids_lt0, pfids_gt0;
     vec<IntDistribution> frag_dist;
     {    String distribs_file = run_dir + "/jump_reads_filt_cpd.distribs";
          vec<IntDistribution> distribs;
          BinaryReader::readFile( distribs_file.c_str( ), &distribs );
          if ( 0 == 1 )
          {
          for ( size_t i = 0; i < distribs.size( ); i++ ) 
          {    IntDistribution dist_lt, dist_gt;
               distribs[i].split( &dist_lt,  &dist_gt );
               dist_lt = dist_lt.reverse( );
               pfids_lt0.push( ProbFuncIntDist(dist_lt) );
               pfids_gt0.push( ProbFuncIntDist(dist_gt) );    }
          }
          String frag_distribs_file = run_dir + "/frag_reads_filt_cpd.distribs";
          BinaryReader::readFile( frag_distribs_file.c_str( ), &frag_dist );    }

     // Go through the contigs.

     vec< triple<int,int,int> > markups;
     vec<Bool> DONE( tigs.size( ), False );
     vec<String> REPORT( tigs.size( ) );
     int REPORT_INDEX = 0;
     uint dots_printed = 0, tigs_processed = 0;
     if ( TIGX < 0 )
     {    cout << Date( ) << ": processing " << tigs.size( ) << " contigs";
          if ( VERBOSITY == 0 ) cout << " (100 dots to follow)" << endl;    }
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int TIG = 0; TIG < (int) tigs.size( ); TIG++ )
     {    if ( TIGX >= 0 && TIG != TIGX ) continue;
          if ( VERBOSITY == 0 && TIGX < 0 )
          {    uint n = tigs.size( );
               uint& i = tigs_processed;
               #pragma omp critical
               {    i++;
                    uint newdots = ((100 * i) / n) - (100 * (i-1)) / n;
                    if (newdots > 0)
                    {    for (uint j = 0; j < newdots; j++)
                         {    cout << ".";
                              dots_printed++;
                              if (dots_printed % 50 == 0) cout << "\n";
                              else if (dots_printed % 10 == 0) cout << " ";    }
                         flush(cout);    }    }    }
          basevector tig = tigs[TIG];
          ostringstream toutx;
          ostream& tout = ( DIRECT ? cout : toutx );
          tout << "\n" << Date( ) << ": analyzing contig " << TIG << " of " 
               << tigs.size( ) << ", l = " << tig.size( ) << " bases" << endl;
          int START = 0, STOP = tig.size( );
          if ( STARTX >= 0 ) START = STARTX;
          if ( STOPX >= 0 ) STOP = STOPX;

          // Load read locations and remove those that are not within bounds.

          vec<read_loc> locs;
          tout << Date( ) << ": loading read locations" << endl;
          #pragma omp critical
          {    locs_file.LoadContig( TIG, locs );    }
          vec<Bool> locs_delete( locs.size( ), False );
          for ( int j = 0; j < locs.isize( ); j++ )
          {    if ( START >= 0 && locs[j].Stop( ) <= START ) locs_delete[j] = True;
               if ( STOP >= 0 && locs[j].Start( ) >= STOP ) 
                    locs_delete[j] = True;    }
          EraseIf( locs, locs_delete );
     
          // Sort out read locations.

          vec<int> ids_frag, ids_jump;
          vec<Bool> or_frag, or_jump;
          vec<int> loc_id_frag, loc_id_jump;
          for ( int i = 0; i < locs.isize( ); i++ )
          {    const read_loc& rl = locs[i];
               if ( rl.Frag( ) ) 
               {    ids_frag.push_back( rl.ReadId( ) );
                    or_frag.push_back( rl.Fw( ) );
                    loc_id_frag.push_back(i);    }
               if ( rl.Jump( ) ) 
               {    ids_jump.push_back( rl.ReadId( ) );
                    or_jump.push_back( rl.Fw( ) );
                    loc_id_jump.push_back(i);    }    }
          SortSync( ids_frag, loc_id_frag );
          SortSync( ids_jump, loc_id_jump );

          // Load reads.

          tout << Date( ) << ": loading reads" << endl;
          vecbasevector bases;
          bases.Read( run_dir + "/frag_reads_filt.fastb", ids_frag );
          bases.Read( run_dir + "/jump_reads_filt.fastb", ids_jump );
          vecqualvector quals;
          quals.Read( run_dir + "/frag_reads_filt.qualb", ids_frag );
          quals.Read( run_dir + "/jump_reads_filt.qualb", ids_jump );
     
          // Put reads in contig orientation.
     
          tout << Date( ) << ": orienting pairs" << endl;
          for ( int i = 0; i < ids_frag.isize( ); i++ )
          {    if ( !locs[ loc_id_frag[i] ].Fw( ) )
               {    bases[i].ReverseComplement( );
                    quals[i].ReverseMe( );    }    }
          for ( int i = 0; i < ids_jump.isize( ); i++ )
          {    if ( !locs[ loc_id_jump[i] ].Fw( ) )
               {    bases[ i + ids_frag.isize( ) ].ReverseComplement( );
                    quals[ i + ids_frag.isize( ) ].ReverseMe( );    }    }
     
          // Find aligned pairs.  Note that this does not yet use innies.

          tout << Date( ) << ": finding aligned pairs" << endl;
          vec< triple<ho_interval,int,int> > fragbounds;
          for ( int i = 0; i < locs.isize( ); i++ )
          {    const read_loc& rl = locs[i];
               if ( !rl.PartnerPlaced() || rl.PartnerContigId() != rl.ContigId() )
                    continue;
               if ( rl.Catywampus( ) || rl.Backwards( ) || rl.Degenerate( ) ) 
                    continue;
               if ( rl.Rc( ) ) continue;
               double offby = Abs( double( rl.SepDelta( ) ) / double( rl.Dev( ) ) );
               if ( offby > heur.max_offby ) continue;
               /*
               if ( rl.Frag( ) ) tout << "(frag) ";
               if ( rl.Jump( ) ) tout << "(jump) ";
               if ( rl.LongJump( ) ) tout << "(long) ";
               tout << rl.Start( ) << "-" << rl.PartnerStop( ) << ", "
                    << "fragment size = " << rl.PartnerStop( ) - rl.Start( ) 
                    << ", off by " << offby << " devs\n";
               */

               /*
               if ( rl.Jump( ) )
               {
               tout << rl.Start( ) << "-" << rl.PartnerStop( ) << ", "
                    << "fragment size = " << rl.PartnerStop( ) - rl.Start( ) 
                    << ", off by " << offby << " devs\n";
               }
               */
     
               // if ( !rl.Jump( ) ) continue; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               fragbounds.push( ho_interval( rl.Start( ), rl.PartnerStop( ) ),
                    rl.SepDelta( ), rl.Dev( ) );    }

          // =======================================================================

          // Troll for indels, using 1 kb window as bait.

          /*
          const int troll_window = 1000;
          const int troll_delta = 100;
          for ( int pos = START; pos < STOP; pos += troll_delta )
          {    int low = pos, high = pos + troll_window;
               int nbounds = 0;
               vec<int> gap, gapdev;
               for ( int j = 0; j < fragbounds.isize( ); j++ )
               {    if ( Subset( ho_interval( low, high ), fragbounds[j].first ) )
                    {    nbounds++;
                         gap.push_back( fragbounds[j].second );
                         gapdev.push_back( fragbounds[j].third );    }    }
               if ( nbounds > 0 )
               {    int gap_ave, gapdev_ave;
                    GapStatsAlt( gap, gapdev, gap_ave, gapdev_ave );
                    PRINT5_TO( tout, low, high, nbounds, 
                         gap_ave, gapdev_ave );    }    }
          */

          // =======================================================================

          // Create data structure to store signal.
          // Current criterion to flag base as bad, one or more of the following:
          // (a) pileup of >= 10 bases at < 50% concordance with reference
          // (b) ambiguity
          // (c) repeat (not implemented).

          tout << Date( ) << ": looking for signal" << endl;
          vec< triple<int,int,String> > signal;
          vec<Bool> weak;

          // Scan to locate ambiguities.

          int pos = 0;
          efasta& T = tigse[TIG];
          for ( int i = 0; i < T.isize( ); i++ )
          {    if ( T[i] == '{' )
               {    int start = pos, stop = pos, istart = i;
                    for ( i++; i < T.isize( ); i++ )
                    {    if ( T[i] == ',' ) break;
                         stop++; pos++;    }
                    for ( i++; i < T.isize( ); i++ )
                         if ( T[i] == '}' ) break;
                    if ( IntervalOverlap( start, stop, START, STOP ) > 0 )
                    {    signal.push( start, stop, "ambiguity "
                              + T.substr( istart, i + 1 - istart ) );    
                         weak.push_back(True);    }    }
               else pos++;    }

          // Scan to locate repeats.
     
          for ( int cpos = 0; cpos < tig.isize( ); cpos++ )
          {    for ( int period = 1; period <= heur.max_tandem_period; period++ )
               {    if ( cpos + heur.min_tandem_copies * period > tig.isize( ) ) 
                         break;
                    int len;
                    for ( len = 1; cpos + len < tig.isize( ); len++ )
                         if ( tig[cpos + len] != tig[cpos + len % period] ) break;
                    if ( len < heur.min_tandem_length ) continue;
                    if ( len < heur.min_tandem_copies * period ) continue;
                    signal.push( cpos, cpos + len, 
                         "tandem repeat, period " + ToString(period) );
                    weak.push_back(True); 
                    cpos += len - 1;    }    }
      
          // Form pileup, identify weak positions and strong windows.

          vec<dumbcall> calls, qcalls;
          calls.resize_and_set( tig.size( ), dumbcall( ) );
          qcalls.resize_and_set( tig.size( ), dumbcall( ) );
          for ( int i = 0; i < ids_frag.isize( ); i++ )
          {    const basevector& b = bases[i];
               const qualvector& q = quals[i];
               AddToPileup( locs[ loc_id_frag[i] ], b, tig, calls );
               AddToPileup( locs[ loc_id_frag[i] ], b, q, tig, qcalls );    }
          for ( int i = 0; i < ids_jump.isize( ); i++ )
          {    const basevector& b = bases[ i + ids_frag.isize( ) ];
               const qualvector& q = quals[ i + ids_frag.isize( ) ];
               AddToPileup( locs[ loc_id_jump[i] ], b, tig, calls );
               AddToPileup( locs[ loc_id_jump[i] ], b, q, tig, qcalls );    }
          for ( int j = START; j < STOP; j++ )
          {    int ncalls = 0, nqcalls = 0, nqcalls_agree = 0;
               for ( int k = 0; k < 6; k++ )
               {    ncalls += calls[j].base[k];
                    nqcalls += qcalls[j].base[k];
                    if ( k == tig[j] ) nqcalls_agree += qcalls[j].base[k];    }
               double qagree = double(nqcalls_agree)/double(nqcalls);
               if ( ncalls >= heur.min_calls && qagree < heur.agree_ceil_weak ) 
               {    signal.push( 
                         j, j+1, "agree = " + ToString( 100.0 * qagree ) + "%" );
                    weak.push_back(True);    }    }
          vec<Bool> marked( STOP - START, False );
          for ( int i = 0; i < signal.isize( ); i++ )
          {    for ( int j = signal[i].first; j < signal[i].second; j++ )
                    if ( j >= START && j < STOP ) marked[ j - START ] = True;    }
          vec<Bool> strong( STOP - START, False );
          for ( int j = START; j < STOP; j++ )
          {    int ncalls = 0, nqcalls = 0, nqcalls_agree = 0;
               for ( int k = 0; k < 6; k++ )
               {    ncalls += calls[j].base[k];
                    nqcalls += qcalls[j].base[k];
                    if ( k == tig[j] ) nqcalls_agree += qcalls[j].base[k];    }
               double qagree = double(nqcalls_agree)/double(nqcalls);
               if ( ncalls >= heur.min_calls && qagree >= heur.agree_floor_strong
                    && !marked[ j - START ] )
               {    strong[ j - START ] = True;    }    }
          for ( int i = 0; i < strong.isize( ); i++ )
          {    if ( !strong[i] ) continue;
               int j;
               for ( j = i + 1; j < strong.isize( ); j++ )
                    if ( !strong[j] ) break;
               if ( j - i >= K ) 
               {    signal.push( START + i, START + j, "" );
                    weak.push_back(False);    }
               i = j - 1;    }

          // Define windows.

          SortSync( signal, weak );
          vec< pair<int,int> > windows;
          vec<String> reports;
          for ( int i = 0; i < signal.isize( ); i++ )
          {    if ( weak[i] )
               {    int j;
                    for ( j = i + 1; j < signal.isize( ); j++ )
                         if ( !weak[j] ) break;
                    if ( i > 0 && j < signal.isize( ) )
                    {    int start = signal[i-1].second - K;
                         int stop = signal[j].first + K;
                         windows.push( start, stop );
                         ostringstream out;
                         for ( int l = i; l < j; l++ )
                         {    out << TIG << "." << signal[l].first << "-" 
                                   << signal[l].second << "   " << signal[l].third 
                                   << "\n";    }
                         reports.push_back( out.str( ) );    }
                    i = j - 1;    }    }

          // Get maximum read length.

          vec<int> readlengths;
          for ( size_t i = 0; i < bases.size( ); i++ )
               readlengths.push_back( bases[i].size( ) );
          UniqueSort(readlengths);
          tout << Date( ) << ": nreadlengths = " << readlengths.size( ) << endl;

          // Go through the windows.

          vec<String> report( windows.size( ) );
          vec< triple<int,int,int> > markupsx( windows.size( ) );
          vec<Bool> markupsx_used( windows.size( ), False );
          vec< triple<int,int,String> > replacements;
          vec< triple<int,int,String> > replacementsx( windows.size( ) );
          vec<Bool> replacementsx_used( windows.size( ), False );
          vec<Bool> done( windows.size( ), False );
          tout << Date( ) << ": traversing " << windows.size( ) << " windows" 
               << endl;
          redo:
          #pragma omp parallel for
          for ( int i = 0; i < windows.isize( ); i++ )
          {    if ( done[i] ) continue;
               if ( WINDOW_START >= 0 && windows[i].first != WINDOW_START ) continue;
               if (DIRECT)
               {    cout << "\n=============================================="
                         << "======================================\n";
                    cout << "\n" << TIG << "." << windows[i].first << "-" 
                         << windows[i].second << "   WINDOW " << i << "\n";    }
               ostringstream routx;
               ostream& rout = ( DIRECT ? cout : routx );
               vec< triple<int,int,int> > markups_this;
               vec< triple<int,int,String> > replacements_this;
               ProcessWindow( heur, logc, head + ".FixLocal." + ToString(i) + "x",
                    windows[i], reports[i], fragbounds, tig, T, TIG, K, ids_frag, 
                    ids_jump, locs, loc_id_frag, loc_id_jump, bases, quals, 
                    readlengths, rout, markups_this, replacements_this );
               if ( !DIRECT ) report[i] = routx.str( );
               if ( markups_this.nonempty( ) ) 
               {    markupsx[i] = markups_this[0];
                    markupsx_used[i] = True;    }
               if ( replacements_this.nonempty( ) )
               {    replacementsx[i] = replacements_this[0];    
                    replacementsx_used[i] = True;    }    }
          tout << Date( ) << ": done traversing windows" << endl;

          // Print reports.

          if ( VERBOSITY >= 2 && !DIRECT )
          {    for ( int i = 0; i < report.isize( ); i++ )
               {    if ( !done[i] ) 
                    {    tout << "\n=============================================="
                              << "======================================\n";
                         tout << "\n" << TIG << "." << windows[i].first << "-" 
                              << windows[i].second << "   WINDOW " << i << "\n";
                         tout << report[i];    }    }    }
          for ( int i = 0; i < report.isize( ); i++ )
               done[i] = True;

          // Check for overlapping windows.  Note that some of this is unnecessary:
          // the original windows overlap in some cases, which didn't have to be
          // the case.  On the other hand, by pushing the boundaries of some 
          // windows, we will occasionally get overlaps that could not have been
          // prevented.

          vec<Bool> to_remove( windows.size( ), False );
          for ( int i = 0; i < windows.isize( ) - 1; i++ )
          {    if ( to_remove[i] ) continue;
               if ( ( markupsx_used[i] || replacementsx_used[i] )
                    && ( markupsx_used[i+1] || replacementsx_used[i+1] )
                    && windows[i].second >= windows[i+1].first )
               {    int low = Min( windows[i].first, windows[i+1].first );
                    int high = Max( windows[i].second, windows[i+1].second );
                    windows[i].first = low, windows[i].second = high;
                    reports[i] = reports[i] + reports[i+1];
                    done[i] = False;
                    to_remove[i+1] = True;    }    }
          int redos = done.isize( ) - Sum(done);
          if ( redos > 0 )
          {    if ( VERBOSITY >= 2 )
               {    tout << "\n" << Date( ) << ": detected " << redos << " overlaps"
                         << " between windows, redoing some" << endl;    }
               EraseIf( windows, to_remove );
               EraseIf( done, to_remove );
               EraseIf( report, to_remove );
               EraseIf( reports, to_remove );
               EraseIf( markupsx, to_remove );
               EraseIf( markupsx_used, to_remove );
               EraseIf( replacementsx, to_remove );
               EraseIf( replacementsx_used, to_remove );
               goto redo;    }

          // Push back results.

          for ( int j = 0; j < windows.isize( ); j++ )
          {    if ( markupsx_used[j] ) markups.push_back( markupsx[j] );
               if ( replacementsx_used[j] ) 
                    replacements.push_back( replacementsx[j] );    }

          // Insert replacements.

          ReverseSort(replacements);
          for ( int i = 0; i < replacements.isize( ); i++ )
          {    int estart = replacements[i].first, estop = replacements[i].second;
               String left = T.substr( 0, estart );
               String middle = replacements[i].third;
               String right = T.substr( estop, T.isize( ) - estop );
               T = left + middle + right;    }

          // Dump report.

          DONE[TIG] = True, REPORT[TIG] = toutx.str( );
          if ( VERBOSITY >= 1 )
          {
               #pragma omp critical
               {    for ( ; REPORT_INDEX < (int) tigs.size( ); REPORT_INDEX++ )
                    {    if ( !DONE[REPORT_INDEX] ) break;
                         cout << REPORT[REPORT_INDEX];    
                         flush(cout);    }    }    }    }

     // Clean up and write new assembly.

     if ( WRITE && TIGX < 0 )
     {    cout << "\n" << Date( ) << ": writing assembly" << endl;
          for ( int s = 0; s < scaffolds.isize( ); s++ )
          {    superb& S = scaffolds[s];
               for ( int j = 0; j < S.Ntigs( ); j++ )
                    S.SetLen( j, tigse[ S.Tig(j) ].Length1( ) );    }
          Assembly A( scaffolds, tigse );
          A.WriteAll( sub_dir + "/" + SCAFFOLDS_OUT );
          Ofstream( out, sub_dir + "/" + SCAFFOLDS_OUT + ".markup" );
          for ( int i = 0; i < markups.isize( ); i++ )
          {    out << markups[i].first << " " << markups[i].second
                    << " " << markups[i].third << "\n";    }    }
     cout << "\n" << Date( ) << ": done, total time used = " << TimeSince(clock) 
          << endl;    }
