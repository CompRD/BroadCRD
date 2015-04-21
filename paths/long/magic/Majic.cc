///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Majic.  Assemble simulated magic reads.
// Files created: sim.fastb, corr.fastb, genome.{fastb,lookup}, 
// genome2.{fastb,lookup}, assembly.{fasta,shbv,dot}.

// revision   genome   edges    end    gaps    error     time           mach
//            size(K)  per kb   gaps           rate(%)
//
// 48583      50       55.5     2      0       1.2       21 seconds     crd5
// 48584      50       12.3     2      0       0.66      4.9 minutes    crd5
// 48589      50       16.8     0      0       0.27      4.0 minutes    crd5
// 48590      50       16.8     0      0       0.27      48 seconds     crd5
// 48591      50       11.2     0      0       0.23      1.65 minutes   crd5
// 48592      50        9.6     0      0       0.19      2.4 minutes    crd5
// 48593      50        7.2     0      0       0.19      2.5 minutes    crd5
// 48594      50        6.8     0      0       0.030     2.4 minutes    crd5
// 48599      50        6.8     0      0       0.030     2.1 minutes    crd8
// 48604      50        6.8     0      0       0.030     1.6 minutes    crd8
// 48607      50        6.6     0      0       0.030     1.6 minutes    crd8
// 48614      50        6.6     0      0       0.030     1.7 minutes    crd5
// 48632      50        6.5     0      0       0.020     1.2 minutes    crd5
// 48637      50        3.3     0      0       0.020     1.2 minutes    crd5

// June 2014 real data test case
// Majic REAL=True TRACE_IDS="[0,100)" TRACE_LEVEL=2 OVER_FOR_TRACE=4 
//      XPASSES=1 MPM=20 DEL1=False

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "PrintAlignment.h"
#include "VecUtilities.h"
#include "efasta/EfastaTools.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "paths/HyperBasevector.h"
#include "paths/KmerBaseBroker.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/EvalAssembly.h"
#include "paths/long/KmerCount.h"
#include "paths/long/LongHyper.h"
#include "paths/long/MakeSimReads.h"
#include "paths/long/Simulation.h"
#include "paths/long/ultra/FounderAlignment.h"
#include "polymorphism/Edit.h"
#include "math/Matrix.h"
#include <omp.h>

struct FriendAlign {
     FriendAlign( int _readid, int _offset, align const& _a, int _errors ) : readid(_readid),
     offset(_offset), a(_a ), errors(_errors) {};

     int readid;
     int offset;         // initial offset
     align a;
     int errors;
};

template<int K1> void FindMatches( const vecbasevector& bases,
     vec< triple<int,int,int> >& match2, vec<int>& pos1s )
{    cout << Date( ) << ": building kmers" << endl;
     vec< triple<kmer<K1>,int,int> > kmers_plus;
     MakeKmerLookup0( bases, kmers_plus );
     cout << Date( ) << ": building matches" << endl;
     int N = bases.size( ) / 2;
     for ( int i = 0; i < kmers_plus.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < kmers_plus.isize( ); j++ )
               if ( kmers_plus[j].first != kmers_plus[i].first ) break;
          for ( int k1 = i; k1 < j; k1++ )
          for ( int k2 = i; k2 < j; k2++ )
          {    if ( k2 == k1 ) continue;
               int id1 = kmers_plus[k1].second, id2 = kmers_plus[k2].second;
               if ( id1 >= N ) continue;
               int pos1 = kmers_plus[k1].third, pos2 = kmers_plus[k2].third;
               int offset = pos1 - pos2;
               match2.push( id1, id2, offset );
               pos1s.push_back(pos1);    }
          i = j - 1;    }
     ParallelSortSync( match2, pos1s );    }

void FindValidated( const vecbasevector& bases,
        int founderid,
        const vec<vec<FriendAlign>>& all_aligns, const int k,
        vecbitvector& validated )
{
     int N = bases.size() / 2;
     const int vs = validated.size();
     ForceAssertEq(vs, N*2);
     ForceAssertLt(founderid, N);

     basevector const& rd2 = bases[founderid];
     validated[founderid].resize(rd2.size(),0);

     for ( auto const& friendalign : all_aligns[founderid] ) {
          basevector const& rd1 = bases[ friendalign.readid ];
          align const& x = friendalign.a;

          int p1 = x.pos1( ), p2 = x.pos2( );
          for ( int j = 0; j < x.Nblocks( ); j++ )
          {    if ( x.Gaps(j) > 0 )
               {    p2 += x.Gaps(j);    }
               else if ( x.Gaps(j) < 0 )
               {    p1 -= x.Gaps(j);    }
               int l = 0;
               while ( l < x.Lengths(j) ) {
                    int lupper;
                    for ( lupper=l ; lupper < x.Lengths(j); lupper++ )
                         if ( rd1[p1+lupper] != rd2[p2+lupper] ) break;
                    if ( lupper - l >= k  ) {
                         for ( ; l < lupper; ++l ) {
                              validated[founderid].Set(p2+l,1);
                         }
                    }
                    ++l;     // we know this is a non-match base or past end
               }
               p1 += x.Lengths(j); p2 += x.Lengths(j);
          }
     }
     validated[founderid+N]=validated[founderid];
     validated[founderid+N].ReverseMe();
}


String Prout( const int id, const int N )
{    if ( id < N ) return ToString(id);
     return ToString( id - N ) + "'";    }

matrix<double>  rates_to_confusion_matrix( double p_d, double p_i, double p_s )
{
     // rates to confusion matrix p(D|T):
     // p_good = 1.0 - p_ind - p_del - p_sub
     // p( base_1 | base_2 ) = p_good      \forall base_1 == base_2
     // p( _ | _ ) = p_good
     // p( base_1 | base_2 ) = p_sub/3.    \forall base_1 != base_2
     // p( base | _ ) = p_ins
     // p( _ | base ) = p_del
     matrix<double> cprob(5,5,0.);

     double p_good = 1.0 - p_i - p_s - p_d;

     double bias = 0.001;       // ensure that nothing has zero probability

     for ( size_t i = 0; i < 5; ++i )
         for ( size_t j = 0; j < 5; ++j ) {
             if ( i == j ) cprob[i][j]=p_good+bias;
             else if ( i == 4 ) cprob[i][j] = p_d+bias;
             else if ( j == 4 ) cprob[i][j] = p_i+bias;
             else cprob[i][j] = p_s / 3.0 + bias;
         }

     // ensure that columns sum to one
     for (size_t j = 0; j < 5; ++j ) {
          double sum = 0.;
          for ( size_t i = 0; i < 5; ++i )
               sum += cprob[i][j];
          for ( size_t i = 0; i < 5; ++i )
               cprob[i][j] /= sum;
     }


     return cprob;
}

void calc_error_rates( VecUCharVec const& multi, UCharVec const& cons,
                        vec<triple<double,double,double>>& rates  )
{
     size_t Nj = multi.size();
     size_t Ni = multi[0].size();

     rates.resize(Nj);
     for ( size_t j = 0; j < Nj; ++j ) {
          // del-ins-sub
          int del = 0;
          int ins = 0;
          int sub = 0;
          double count = 0.;
          for ( size_t i = 0; i < Ni; ++i ) {
               auto call = multi[j][i];
               if ( call != 5 && cons[i] != 5 ) {
                    count+=1.;
                    if ( call < 4 && cons[i] < 4 && call != cons[i] )
                         sub++;
                    else if ( call < 4 && cons[i] == 4 )
                         ins++;
                    else if ( cons[i] < 4 && call == 4 )
                         del++;
               }
          }

          if ( count > 0 ) rates[j] = make_triple( del / count, ins / count, sub / count );
          else rates[j] = make_triple(0.,0.,0.);
     }
}

void bayes_consensus(vec<matrix<double>> const& p_d_t, VecUCharVec const& multi,
          UCharVec& cons)
// note cons is both input and output -- it's expected that an esimtated
// consensus is here simply for calculation of the prior
{
     // consensus:
     //     Pr{ T_i = t | D; \theta } = p(t) \prod_j p(D_j | T=t )
     // ----------------------------------------------------------
     // \sum_t' Pr{ T_i = t' | D; \theta } = p(t') \prod_j p(D_j | T=t' )
     size_t Nj = multi.size();
     size_t Ni = multi[0].size();

     vec<double> priors(5,0.);
     // global prior based on the current estimated consensus
     for ( size_t i = 0; i < Ni; ++i )
          if (cons[i] < 5) priors[cons[i]] += 1.;
     double sum = Sum(priors);
     ForceAssertGe( sum, 1.0 );
     for ( auto& prior : priors ) prior = log(prior/sum);


     for ( size_t i = 0; i < Ni; ++i ) {

          vec<double> calls = priors;
          for ( size_t j = 0; j < Nj; ++j ) {
               for ( size_t k = 0; k < 5; ++k )
                    if ( multi[j][i] < 5 )
                         calls[k] += log( p_d_t[j][ multi[j][i] ][k] );
          }
          vec<int> ident(calls.size(), vec<int>::IDENTITY);
          if ( i == 0 ) {
              /*
              cout << "COL0: ";
              for ( size_t k = 0; k < calls.size(); ++k )
                  cout << exp(calls[k]) <<  " ";
              cout << endl;
              */
          }
          ReverseSortSync(calls, ident );
          if ( i == 0 ) {
              /*
              cout << "sorted COL0: ";
              for ( size_t k = 0; k < calls.size(); ++k )
                  cout << exp(calls[k]) <<  " ";
              cout << endl;
              cout << "sorted log COL0: ";
              for ( size_t k = 0; k < calls.size(); ++k )
                  cout << calls[k] <<  " ";
              cout << endl;
              cout << "sorted ident COL0: ";
              for ( size_t k = 0; k < calls.size(); ++k )
                  cout << ident[k] <<  " ";
              cout << endl;
              */
          }
          cons[i] = ident[0];
     }
}

// begs to be in a class...
typedef enum { CONSENSUS_NONE, CONSENSUS_MAJORITY,
     CONSENSUS_STAT, CONSENSUS_ITER } ConsensusType;

ConsensusType parseConsensusType( String const& s )
{
     if ( s == "MAJORITY" ) return CONSENSUS_MAJORITY;
     else if ( s == "STAT" ) return CONSENSUS_STAT;
     else if ( s == "ITER" ) return CONSENSUS_ITER;
     else if ( s != "NONE" )
          FatalErr("unrecognized consensus type: " + s);

     return CONSENSUS_NONE;

}

String consensusTypeStr( ConsensusType const& ctype )  {
     if ( ctype == CONSENSUS_MAJORITY) return "MAJORITY";
     else if ( ctype == CONSENSUS_STAT ) return "STAT";
     else if ( ctype == CONSENSUS_ITER ) return "ITER";
     else if ( ctype != CONSENSUS_NONE)
          FatalErr("unrecognized consensus code: " + ctype);

     return "NONE";
}

void Consensus( double p_d, double p_i, double p_s, VecUCharVec const& multi,
          size_t founder, vec<triple<double,double,double>>& rates, UCharVec& cons,
           ConsensusType mode = CONSENSUS_MAJORITY, vec<double>* pWeights = 0)
{
     if ( multi.size() == 0 ) {
          cons.clear();
          rates.clear();
          return;
     }
     ForceAssertLt(founder, multi.size() );

     size_t Nj = multi.size();
     size_t Ni = multi[0].size();
     rates.resize( Nj );
     cons.resize( Ni );

     // initialization value for consensus
     // majority vote
     for ( size_t i = 0; i < Ni; ++i )  {
          vec<double> score(6,0.);
          vec<int> ident(6, vec<int>::IDENTITY);
          for ( size_t j = 0; j < Nj; ++j ) {
               ForceAssertLt( multi[j][i], 6 );
               if ( pWeights == 0 ) score[ multi[j][i] ] += 1.0;
               else score[ multi[j][i] ] += (*pWeights)[j];
          }
          ReverseSortSync(score, ident );
          // don't score missing bases, unless they're all missing bases
          // 5 == missing, 4 == gap, <4 == base
          double top_score = score[0];
          int top_base = ident[0];
          size_t top_idx = 0;
          if ( top_base == 5 && score[1] != 0 ) {
               top_score = score[1];
               top_base = ident[1];
               top_idx = 1;
          }
          cons[i] = top_base;           // set provisionally to first high scorer...
          // ... but we make any ties go to the founder base
          for ( size_t k = top_idx; k < score.size() && score[k] == top_score; ++k )
               if ( ident[k] == multi[founder][i] )
                    cons[i] = ident[k];
     }

     // initialization values for rates
     calc_error_rates( multi, cons, rates );

     if ( mode == CONSENSUS_MAJORITY ) return;  // EARLY RETURN

     for ( size_t iter = 0; iter < 10; ++iter ) {
     vec<matrix<double>> p_d_t( Nj );
     for ( size_t j = 0; j < Nj; ++j ) {
         p_d_t[j] = rates_to_confusion_matrix(rates[j].first, rates[j].second, rates[j].third);
         // debugging below
         /*
         PRINT4(j, rates[j].first, rates[j].second, rates[j].third );
         cout << "matrix " << j << ":" << endl;
         for ( int k1 = 0; k1 < p_d_t[j].Nrows(); ++k1 ) {
             for ( int k2 = 0; k2 < p_d_t[j].Ncols(); ++k2 )
                 cout << p_d_t[j][k1][k2] << " ";
             cout << endl;
         }
         */

     }

     bayes_consensus( p_d_t, multi, cons );

     if ( mode == CONSENSUS_STAT ) return; //EARLY RETURN
     calc_error_rates( multi, cons, rates );
     }
}

void write_uchar_reads( ostream& out, String const& label, UCharVec read )
{
    size_t count = 0;
    out << ">" << label << endl;
    for (size_t i = 0; i < read.size(); ++i ) {
        if ( read[i] < 4 ) {
            out << as_base( read[i] );
            if ( ++count % 80 == 0 ) out << endl;
        }
    }
    if ( count % 80 != 0 ) out << endl;
}

void MakeAlignments( const basevector& b, const vecbasevector& bases,
     const vec< triple<int,int,int> >& all_over, const int min_votes,
     const int bandwidth, const int mpm, ostream& out, const int trace_level,
     const int id1, const int over_for_trace, vec<vec<FriendAlign>>& all_aligns )
{

     int N = bases.size( ) / 2;
     vec< vec<edit0> > edits( b.size( ) );
     for ( int i = 0; i < all_over.isize( ); i++ )
     {    align x;
          int errors;
     
          SmithWatBandedA( bases[ all_over[i].first ], b,
               -all_over[i].second, bandwidth, x, errors, 0, 1, 1 );    
     
          if ( trace_level >= 1 )
          {    out << "\nalignment of read " << id1 << " to read ";
               if ( all_over[i].first < N ) out << all_over[i].first;
               else out << all_over[i].first - N << "'";
               out << ", score=" << errors << endl;
               PrintVisualAlignment( True, out, bases[ all_over[i].first ], b, x );
               vec<ho_interval> perf1;
               x.PerfectIntervals1( bases[ all_over[i].first ], b, perf1 );
               int M = 0;
               for ( int l = 0; l < perf1.isize( ); l++ )
                    M = Max( M, perf1[l].Length( ) );
               out << "max perfect match = " << M << "\n";
               if ( M < 16 ) 
                    out << "Warning: seems like a crappy alignment." << endl;     }

          all_aligns[id1].emplace_back( all_over[i].first, all_over[i].second, x, errors );
     }
     
}

basevector squash_gaps( UCharVec const& in )
{
     ostringstream s;
     for ( auto const& base : in  )
          if ( base < 4 ) s << as_base(base);

     return basevector(s.str());
}

basevector InferEdits( const basevector& b, const vecbasevector& bases,
     vecString& bases_names, const vec<vec<FriendAlign>>& all_aligns, 
     vec<vec<double>>& all_align_scores, const int min_votes, const int mpm, 
     ostream& out, const int trace_level, const int id1, const int over_for_trace, 
     const Bool del1, double max_friend_score, ConsensusType ctype )
{
     // Generate multiple alignment.  Note redundant generation of pairwise 
     // alignments.

     if ( ctype != CONSENSUS_NONE || trace_level >= 2 )
     {    const double del_rate = 0.05;
          const double ins_rate = 0.02;
          const double sub_rate = 0.05;
          Scorer scorer( sub_rate, del_rate, ins_rate );
          VecUCharVec multi;
          vecbasevector gang;
          vec<double> scores;
          gang.push_back(b);
          vec<align> aligns;
          aligns.push_back( align( ) );
          scores.push_back(0.);         // dummy score for the founder read
          vec<int> ids;
          ids.push_back(id1);
          for ( int i = 0; i < all_aligns[id1].isize( ); i++ )
          {    int friendid = all_aligns[id1][i].readid;
               align x = all_aligns[id1][i].a;
               int errors = all_aligns[id1][i].errors;
               // filter false friends based on max_friend_score
               double friend_score = all_align_scores[id1][i];
               if ( friend_score <= max_friend_score)  {
                    vec<ho_interval> perf1;
                    x.PerfectIntervals1( bases[ friendid ], b, perf1 );
                    int M = 0;
                    for ( int l = 0; l < perf1.isize( ); l++ )
                         M = Max( M, perf1[l].Length( ) );
                    if ( M < mpm ) continue;
                    x.Flip( );
                    gang.push_back( bases[ friendid ] );
                    aligns.push_back(x);
                    ids.push_back( friendid );
                    scores.push_back( friend_score );
               }
          }
          if ( ctype != CONSENSUS_NONE || aligns.isize( ) - 1 >= over_for_trace )
          {

          AlignFriendsToFounder( gang, 0, aligns, scorer, &multi );

          out << "\nmultiple alignment of read " << id1 << "\n\n";
          out << "to reads";
          for ( int i = 0; i < ids.isize( ); i++ )
               out << " " << Prout( ids[i], bases.size( )/2 );
          out << "\n\n";

          // Delete columns having only one base, and at least five entries.

          if (del1)
          {    vec<Bool> to_delete( multi[0].size( ), False );
               for ( int j = 0; j < (int) multi[0].size( ); j++ )
               {    int non45 = 0;
                    for ( int i = 0; i < (int) multi.size( ); i++ )
                         if ( multi[i][j] != 4 && multi[i][j] != 5 ) non45++;
                    if ( non45 == 1 && multi.size( ) >= 5 ) to_delete[j] = True;    }
               for ( int i = 0; i < (int) multi.size( ); i++ )
               {    SerfVec<uchar> m = multi[i];
                    vec<uchar> mx;
                    for ( int l = 0; l < (int) m.size( ); l++ )
                    mx.push_back( m[l] );
                    EraseIf( mx, to_delete );
                    multi[i].resize(0);
                    for ( int l = 0; l < mx.isize( ); l++ )
                         multi[i].push_back( mx[l] );    }    }

          // Print the multi-alignment.

          UCharVec cons1,cons2;
          vec<triple<double,double,double>> rates;

          // weights are for majority voting (CONSENSUS_MAJORITY) and
          // for initialization of other methods (which are initialized
          // with majority vote).  A weight of 2, will score as if the
          // read appeared twice in the stack.

          vec<double> weights( multi.size(), 1.0);
          Consensus( del_rate, ins_rate, sub_rate, multi, 0, rates, cons1, 
               CONSENSUS_MAJORITY, &weights );
          Consensus( del_rate, ins_rate, sub_rate, multi, 0, rates, cons2, 
               ctype );
          multi.push_back( cons1 );
          multi.push_back( cons2 );

#if 0
          Ofstream( debug_out, "neil.fasta" );
          write_uchar_reads(debug_out, "founder", multi[0] );
          write_uchar_reads(debug_out, "majority", cons1);
          write_uchar_reads(debug_out, "stat", cons2);
#endif


          if ( trace_level >= 2 )
          {
          const int width = 80;
          for ( int js = 0; js < (int) multi[0].size( ); js += width )
          {    for ( int i = 0; i < (int) multi.size( ); i++ )
               {    Bool nothing = True;
                    for ( int j = js; j < Min( (int) multi[0].size( ), js + width ); j++ )
                    {    uchar c = multi[i][j];
                         if ( c != 5) nothing = False; }
                    if ( !nothing )
                    {
                         for ( int j = js; j < Min( (int) multi[0].size( ), js + width ); j++ )
                         {    uchar c = multi[i][j];
                              if ( c < 4 ) out << as_base(c);
                              if ( c == 4 ) out << "-";
                              if ( c == 5 ) out << "=";  }

                         if ( i < (int)multi.size()-2 ) {

                              if ( bases_names.size( ) > 0 )
                                   out << " " << bases_names[ ids[i] ];

                              triple<double,double,double> 
                                   rate=make_triple(-1,-1,-1);
                              rate = rates[i];
                              out << " [" << (int)(rate.first*100) << ","
                              << (int)(rate.second*100)
                              << "," << (int)(rate.third*100) << "] (" << scores[i] << ")";
                         } else if ( i == (int)multi.size()-2 ) out << " [majority]" ;
                         else if ( i == (int)multi.size()-1 ) out << " [" + consensusTypeStr(ctype) << "]" ;
                         out << "\n";  }  }
               out << "\n";    }
          }
          if ( ctype != CONSENSUS_NONE ) return squash_gaps( cons2 );    // EARLY RETURN!
          } }

     // Proceed.

     int N = bases.size( ) / 2;
     vec< vec<edit0> > edits( b.size( ) );
     for ( int i = 0; i < all_aligns[id1].isize( ); i++ )
     {
          int friendid = all_aligns[id1][i].readid;
          align x = all_aligns[id1][i].a;
          int errors = all_aligns[id1][i].errors;

          if ( trace_level >= 1 )
          {    out << "\nalignment of read " << id1 << " to read ";
               if ( friendid <  N ) out << friendid;
               else out << friendid - N << "'";
               out << ", score=" << errors << endl;
               PrintVisualAlignment( True, out, bases[ friendid ], b, x );
               vec<ho_interval> perf1;
               x.PerfectIntervals1( bases[ friendid ], b, perf1 );
               int M = 0;
               for ( int l = 0; l < perf1.isize( ); l++ )
                    M = Max( M, perf1[l].Length( ) );
               out << "max perfect match = " << M << "\n";
               if ( M < 16 ) 
                    out << "Warning: seems like a crappy alignment." << endl;     }
     
          basevector rd1 = bases[ friendid ];
          basevector rd2 = b;
     
          int p1 = x.pos1( ), p2 = x.pos2( );
          for ( int j = 0; j < x.Nblocks( ); j++ )
          {    if ( x.Gaps(j) > 0 ) 
               {    out << "EDIT: delete " 
                         << x.Gaps(j) << " bases at " << p2 << endl;
                    edits[p2].push( DELETION, x.Gaps(j) );
                    p2 += x.Gaps(j);    }
               if ( x.Gaps(j) < 0 ) 
               {    out << "EDIT: insert ";
                    String s;
                    for ( int l = 0; l < -x.Gaps(j); l++ )
                    {    out << as_base( rd1[ p1 + l ] );
                         s.push_back( as_base( rd1[ p1 + l ] ) );    }
                    out << " at " << p2 << endl;
                    edits[p2].push( INSERTION, s );
                    p1 -= x.Gaps(j);    }
               for ( int l = 0; l < x.Lengths(j); l++ )
               {    if ( rd1[p1] != rd2[p2] )
                    {    out << "EDIT: change to " << as_base( rd1[p1] )
                              << " at " << p2 << endl;
                         edits[p2].push( SUBSTITUTION, as_base( rd1[p1] ) );    }
                    ++p1; ++p2;    }    }    }
     
     for ( int l = 0; l < b.isize( ); l++ )
          Sort( edits[l] );
     /*
     out << "\nSORTED:\n\n";
     for ( int l = 0; l < b.isize( ); l++ )
     {    for ( int m = 0; m < edits[l].isize( ); m++ )
          {    int n = edits[l].NextDiff(m);
               if ( n - m >= 2 )
                    out << l << " " << edits[l][m] << " [" << n-m << "]" << "\n";
               m = n - 1;    }    }
     */
     basevector mod( b );

     for ( int l = b.isize( ) - 1; l >= 0; l-- )
     {    
          vec< pair<int,edit0> > votes;

          for ( int m = 0; m < edits[l].isize( ); m++ )
          {    int n = edits[l].NextDiff(m);
               votes.push( n - m, edits[l][m] );
               m = n - 1;    }
          ReverseSort(votes);

          for ( int m = 0; m < votes.isize( ); m++ )
          {    if ( votes[m].first >= min_votes )

               {    const edit0& e = votes[m].second;

                    if ( e.etype == SUBSTITUTION ) mod.Set( l, as_char( e.seq[0] ) );

                    if ( e.etype == DELETION )
                    {    basevector mod2;
                         for ( int z = 0; z < mod.isize( ); z++ )
                         {    if ( z >= l && z < l + e.n ) continue;
                              mod2.push_back( mod[z] );    }
                         mod = mod2;    }

                    if ( e.etype == INSERTION )
                    {    basevector mod2;
                         for ( int z = 0; z < mod.isize( ); z++ )
                         {    mod2.push_back( mod[z] );
                              if ( z == l - 1 )
                              {    for ( int w = 0; w < e.seq.isize( ); w++ )
                                   {    mod2.push_back( as_char( 
                                             e.seq[w] ) );    }    }    }
                         mod = mod2;    }

                    break;
     
                    // out << l << " " << e << " [" << n-m << "]" << "\n";
                         }    }    }

     return mod;    
}


void ScoreAlignments( vecbasevector const& bases, vecbitvector const& validated,
          vec<vec<FriendAlign>> const& all_aligns, vec<vec<double>>& all_align_scores )
{
     ForceAssertEq(bases.size()/2,all_aligns.size() );
     const int gap_open_penalty = 12;
     const int gap_extend_penalty = 1;
     const int mismatch_penalty = 3;
     const int min_valid = 2;           // how many bases into valid overlap to count an event
     const int min_scorable = 100;      // don't score a read with fewer than this many valid overlapping positions

     all_align_scores.resize( all_aligns.size() );

#pragma omp parallel for
     for (size_t id = 0; id < bases.size() / 2; ++id ) {

          all_align_scores[id].clear();

          basevector const& rd2 = bases[ id ];

          for ( auto const& al : all_aligns[id] ) {

               auto const& x = al.a;
               basevector const& rd1 = bases[ al.readid ];

               int nvalid = 0;
               int score = 0;
               int nscorable = 0;
               int p1 = x.pos1( ), p2 = x.pos2( );

               for ( int j = 0; j < x.Nblocks( ); j++ ) {
                    if ( x.Gaps(j) > 0 ) {

                         if ( validated[al.readid][p1] ) { nvalid++; }
                         else nvalid = 0;
                         for ( int l = 0; l < x.Gaps(j); ++l ) {
                              if ( !validated[id][p2] ) nvalid = 0;
                              p2++;
                         }
                         if ( nvalid > min_valid ) {
                              nscorable++;
                              score += gap_open_penalty + (x.Gaps(j)-1)*gap_extend_penalty;
                         }
                    }
                    if ( x.Gaps(j) < 0 ) {

                         if ( validated[id][p2] ) { nvalid++; }
                         else nvalid = 0;
                         for ( int l = 0; l < -x.Gaps(j); ++l ) {
                              if ( !validated[al.readid][p1] ) nvalid = 0;
                              p1++;
                         }
                         if ( nvalid > min_valid ) {
                              nscorable++;
                              score += gap_open_penalty + (-x.Gaps(j)-1)*gap_extend_penalty;
                         }

                    }
                    for ( int l = 0; l < x.Lengths(j); l++ ) {
                         if ( validated[id][p2] && validated[al.readid][p1] ) { nvalid++; nscorable++; }
                         else nvalid = 0;
                         if ( nvalid > min_valid ) {
                              nscorable++;
                              if ( rd1[p1] != rd2[p2] ) {
                                   score += mismatch_penalty;
                              }
                          }
                         ++p1; ++p2;
                    }

               } //for each block

               if ( nscorable > min_scorable ) 
                    all_align_scores[id].push_back(score/static_cast<double>(nscorable));
               else
                    all_align_scores[id].push_back(std::numeric_limits<double>::max());
          } // for each align

     } // for each id1

}

int main(int argc, char *argv[])
{
     RunTime( );
     double clock = WallClockTime( );

     BeginCommandArgumentsAcceptEmptyArgListNoHeader;
     CommandArgument_Int_OrDefault_Doc(XPASSES, 4, 
          "number of error correction passes");
     CommandArgument_Int_OrDefault_Doc(START, 0, "start of region");
     CommandArgument_Int_OrDefault_Doc(NB, 50000, "number of bases in region");
     CommandArgument_Bool_OrDefault_Doc(DELETE_LOW_COV, True, 
          "delete low-coverage edges");
     CommandArgument_String_OrDefault_Doc(TRACE_IDS, "", "provide info on these read ids");
     CommandArgument_Int_OrDefault_Doc(TRACE_LEVEL, 1, "1 or 2, for TRACE_ID");
     CommandArgument_Int_OrDefault_Doc(OVER_FOR_TRACE, 0, 
          "minimum number of overlaps to trace a read");
     CommandArgument_Bool_OrDefault_Doc(REAL, False, "use real data instead");
     CommandArgument_Int_OrDefault_Doc(MPM, 0, 
          "minimum perfect match for read-read alignments");
     CommandArgument_Bool_OrDefault_Doc(DEL1, True, 
          "delete columns in multiple alignment if only one base call is present");
     CommandArgument_String_OrDefault_Doc(CONSENSUS, "NONE", 
          "try stack-based consensus: NONE (default), MAJORITY, STAT, "
          "ITER (very experimental)");
     CommandArgument_Double_OrDefault_Doc(MAX_FRIEND_SCORE, std::numeric_limits<double>::max(),
          "use a value (recommend 1.0) to filter friends with a score greater than this");
     CommandArgument_Bool_OrDefault(USE_OLD_LRP_METHOD,True);
     EndCommandArguments;
     uint NUM_THREADS = 48;
     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads(NUM_THREADS);

     // Parse TRACE_IDS.

     vec<int> trace_ids;
     ParseIntSet( TRACE_IDS, trace_ids );

     // Define reference, a circularized chunk of E. coli.

     basevector refb;
     {    vecbasevector G( "/wga/dev/references/Escherichia_coli/genome.fastb" );
          refb = basevector( G[0], START, NB );    }
     vecbasevector genome, genome2;
     genome.push_back(refb);
     genome.WriteAll( "genome.fastb" );
     genome2.push_back( Cat( refb, refb ) );
     genome2.WriteAll( "genome2.fastb" );
     SystemSucceed( "MakeLookupTable SOURCE=genome.fastb OUT_HEAD=genome "
          "LO=True NH=True" );
     SystemSucceed( "MakeLookupTable SOURCE=genome2.fastb OUT_HEAD=genome2 "
          "LO=True NH=True" );

     // Create simulated data.  It is 24x coverage by ~5 kb reads with error
     // profile given by 5% deletions, 2% insertions, 5% substitutions.

     vecbasevector bases;
     vecString bases_names;
     if ( !REAL )
     {    basevector GG;
          for ( int j = 0; j < 10; j++ )
               GG = Cat( GG, refb );
          ref_data ref;
          ref.G.push_back( GG );
          ref.G3 = ref.G;
          ref.is_circular.resize( 1, False );
          ref.ploidy.resize( 1, 1 );
          CreateGlocs( ref.G, ref.LG, ref.Glocs );
          CreateGlocs( ref.G3, ref.LG, ref.G3locs );
          CreateGlocs( ref.G3plus, ref.LG, ref.G3pluslocs );
          String SIM 
               = "{ERR_DEL=0.05,ERR_INS=0.02,ERR_SUB=0.05,COVERAGE=2.4,LEN=5K}";
          vec<ref_loc> readlocs;
          long_sim sim(SIM);
          MakeSimReads( sim, ref, bases, readlocs );    }
     else 
     {    
          FetchReads( bases, bases_names,
               "/wga/scr4/vendor/on/2014-07-09/ecoli_r6/ecoli_2d.select1.fasta" );
          for ( int i = 0; i < (int) bases_names.size( ); i++ )
               bases_names[i] = bases_names[i].Before( " " );

          /*
          vecbasevector bases1, bases2;
          vecString bases_names1, bases_names2;
          FetchReads( bases1, bases_names1, 
               "/wga/scr4/macro/scardo/scardo.Feb7.1.subset1.fasta" );
          FetchReads( bases2, bases_names2, 
               "/wga/scr4/macro/scardo/scardo.Feb11.1.subset1.fasta" );
          bases = bases1;
          bases.Append(bases2);
          bases_names = bases_names1;
          bases_names.Append(bases_names2);    
          */

          }

     bases.WriteAll( "sim.fastb" );
     int N = bases.size( );

     // Heuristics.

     const int K1 = 16;
     const int K2 = 400;
     const int min_spread = 50;
     const int max_offset_diff = 100;
     const int min_start = 50;
     const int min_predicted_overlap = 1000;
     const int min_votes = 10;
     const int min_for_valid = 12; // validated base is in perfect match of this size

     // Carry out XPASSES passes.

     for ( int xpass = 1; xpass <= XPASSES; xpass++ )
     {    
          // Find K1-mer matches between the reads.

          vecbasevector bases_new(bases);
          vecbasevector basesrc(bases);
          for ( int i = 0; i < (int) basesrc.size( ); i++ )
               basesrc[i].ReverseComplement( );
          bases.Append(basesrc);

          vec< triple<int,int,int> > match2;
          vec<int> pos1s;
          if ( xpass == 1 ) FindMatches<K1>( bases, match2, pos1s );
          else FindMatches<100>( bases, match2, pos1s );

          // Use the matches to define overlaps between the reads.  First for the
          // kmer matches between two reads, we define the spread.  We consider
          // all pairs of matches.  If their offsets differ by more than 
          // max_offset_diff, we ignore the pair.  If the start position on the
          // first read is less than min_start, we ignore the pair.  Otherwise the 
          // spread is the maximum over all pairs of the distance between the start
          // positions on the first read.  If the spread is less than min_spread,
          // we ignore the overlaps between the two reads.
          //
          // Otherwise we define an overlap by taking the median offset.

          cout << Date( ) << ": start main loop" << endl;
          vec< vec< triple<int,int,int> > > over(N);
          int overlaps = 0;

          // this is needed to keep from perturbing results, so we still can pick
          // the middle offset, sorted within a range i..j (in match2)
          // match here is sorted by pos1s BEFORE offset

          vec< quad<int,int,int,int> > match;
          for ( size_t i = 0; i < match2.size(); ++i )
               match.push( match2[i].first, match2[i].second, pos1s[i], match2[i].third );
          Sort(match);

          int i;
          for ( i = 0; i < match.isize(); ++i ) {

               // define j so that [i,j) are where id1 and id2 are the same as match[i]
               int id1 = match[i].first, id2 = match[i].second, j;
               for ( j = i + 1; j < match.isize( ); j++ )
                    if ( match[j].first != id1 || match[j].second != id2 ) break;

               int spread = 0;

               // adjust k1 upwards so that pos1 is at least min_start

               int k1;
               for ( k1 = i; k1 < j; k1++ ) if ( match[k1].third >=  min_start ) break;
               if ( k1 != j ) {
                    // run through k1 from bottom
		   for ( ; k1 < j; k1++ ) {
		      // find highest viable k2 at the top
		      // (sorted, so it will score highest)
		      int k2;
		      for ( k2 = j-1; k2 > k1; k2-- )
			   if ( Abs( match[k2].fourth - match[k1].fourth )
				     <= max_offset_diff ) break;
		      if ( k2 == k1 ) continue;
		      spread = Max( spread, Abs( match[k1].third - match[k2].third ) );
		   }
               }


               if ( spread >= min_spread )  {
		    int offset = match2[ i + (j-i)/2 ].third;
		    over[id1].push( id2, offset, j-i );

		    overlaps++;

               }
               i = j - 1;    }

          cout << "\n";
          cout << Date( ) << ": done with main loop" << endl;

          // Generate transitive overlaps.

          vec<vec<FriendAlign>> all_aligns(N);
          vecbitvector validated(2*N);          // we rc the validation info, too
          vec<basevector> corrected;
          vec<String> reports(N);
          #pragma omp parallel for
          for ( int id1 = 0; id1 < N; id1++ )
          {    vec< triple<int,int,int> > all_over = over[id1];
               ostringstream out;
               out << "\noverlaps with id1 = " << id1
                    << ", len = " << bases[id1].size( ) << "\n";
               out << "direct:\n";
               vec<int> ids;
               ids.push_back(id1);
               for ( int i = 0; i < over[id1].isize( ); i++ )
               {    int id2 = over[id1][i].first;
                    ids.push_back(id2);
                    out << "id2 = ";
                    if ( id2 < N ) out << id2;
                    else out << id2 - N << "'";
                    out << ", len = " << bases[id2].size( )
                         << ", offset = " << over[id1][i].second
                         << ", multiplicity = " << over[id1][i].third << "\n";    }
               Sort(ids);
               out << "transitive:\n";
               vec< pair<int,int> > trans;
               for ( int i = 0; i < over[id1].isize( ); i++ )
               {    int id2 = over[id1][i].first;
                    int offset2 = over[id1][i].second;
                    Bool rc2 = False;
                    if ( id2 >= N )
                    {    id2 -= N;
                         rc2 = True;    }
                    for ( int j = 0; j < over[id2].isize( ); j++ )
                    {    int id3 = over[id2][j].first;
                         int offset3 = over[id2][j].second;
                         if (rc2)
                         {    if ( id3 < N ) id3 += N;
                              else id3 -= N;    
                              offset3 = -offset3 
                                   - bases[id3].isize( ) + bases[id2].isize( );    }
                         trans.push( id3, offset2 + offset3 );    
                         if ( BinMember( trace_ids, id1 ) && !BinMember( ids, id3 ) ) 
                         {    out << "xpass " << xpass << ", adding read "
                                   << Prout(id3,N) << " by transitivity with " 
                                   << Prout(over[id1][i].first,N) 
                                   << endl;    }    }    }
               Sort(trans);
               for ( int t = 0; t < trans.isize( ); t++ )
               {    int u;
                    for ( u = t + 1; u < trans.isize( ); u++ )
                         if ( trans[u].first != trans[t].first ) break;
                    int offset = 0;
                    for ( int l = t; l < u; l++ )
                         offset += trans[l].second;
                    offset /= (u-t);
                    int id3 = trans[t].first;
                    if ( !BinMember( ids, id3 ) )
                    {    int pred = IntervalOverlap( 0, bases[id1].isize( ),
                              offset, offset + bases[id3].isize( ) );
                         if ( pred >= min_predicted_overlap )
                         {    out << "id3 = " << Prout(id3,N) << ", len = " 
                                   << bases[id3].size( ) << ", offset = " << offset 
				   << ", multiplicity = " << u - t
                                   << ", predicted overlap = " << pred << "\n";
                              if ( BinMember( trace_ids, id1 ) )
                              {    out << "xpass " << xpass << ", adding read "
                                        << Prout(id3,N) << " with offset " 
                                        << offset << endl;    }
                              all_over.push( id3, offset, u-t );    }    }
                    t = u - 1;    }

               // Make gang and alignments.

               const int bandwidth1 = 200;
               int trace_level = 0;
               if ( BinMember( trace_ids, id1 ) 
                    && all_over.isize( ) >= OVER_FOR_TRACE )
               {    trace_level = TRACE_LEVEL;    }
               if ( trace_level >= 1 )
               {    out << "xpass " << xpass << ", alignments of read "
                         << id1 << endl;
                    for ( int i = 0; i < all_over.isize( ); i++ )
                    {    int id2 = all_over[i].first, offset = all_over[i].second;
                         out << "id2 = ";
                         if ( id2 < N ) out << id2;
                         else out << id2 - N << "'";
                         out << ", offset = " << offset << endl;    }    }

               MakeAlignments( bases[id1], bases, all_over, min_votes, bandwidth1, MPM, out,
                       trace_level, id1, OVER_FOR_TRACE, all_aligns );
               // validated is a bitvector for each read (one bit per
               // base) with "true" if the base participated in a
               // perfect alignment of length "min_for_valid" with at
               // least one other read.
               FindValidated( bases, id1, all_aligns, min_for_valid, validated );

               if (trace_level > 1 ) {
#pragma omp critical
                   cout << out.str( ) << endl; } }

          // Okay, now we have all versus all alignments and the
          // "validated" bitvec for each read.  So let's generate a
          // quality measure for each pairwise alignment scored only
          // using validated bases.  The hope here is that false friends
          // will have wildly different scores.
          vec<vec<double>> all_align_scores;
          ScoreAlignments( bases, validated, all_aligns, all_align_scores );


          #pragma omp parallel for
          for ( int id1 = 0; id1 < N; id1++ )
          {
               ostringstream out;
               int trace_level = 0;
               if ( BinMember( trace_ids, id1 )
                    && all_aligns[id1].isize( ) >= OVER_FOR_TRACE )
               {    trace_level = TRACE_LEVEL;    }
               auto consensus = parseConsensusType( CONSENSUS );
               basevector mod1 = InferEdits( bases[id1], bases, bases_names, 
                    all_aligns, all_align_scores, min_votes, MPM, out, trace_level, 
                    id1, OVER_FOR_TRACE, DEL1, MAX_FRIEND_SCORE, consensus );

               bases_new[id1] = mod1;
               if (trace_level > 1 ) {
#pragma omp critical
                    cout << out.str( ) << endl;
               }
          }


          bases = bases_new;    }
     cout << Date( ) << ": done correcting reads" << endl;
     bases.WriteAll( "corr.fastb" );

     // Form initial assembly.

     VecEFasta correctede;
     for ( int i = 0; i < (int) bases.size( ); i++ )
          correctede.push_back( efasta( bases[i] ) );
     vec<pairing_info> cpartner;
     SupportedHyperBasevector shb;
     long_logging logc( "", "" );
     ref_data ref;
     ref.G.push_back( refb );
     ref.G3 = ref.G;
     ref.is_circular.resize( 1, True );
     ref.ploidy.resize( 1, 1 );
     CreateGlocs( ref.G, ref.LG, ref.Glocs );
     CreateGlocs( ref.G3, ref.LG, ref.G3locs );
     CreateGlocs( ref.G3plus, ref.LG, ref.G3pluslocs );
     vec<basevector> p;
     p.push_back( genome[0] );
     const int KX = 80; // no justification for this
     ref.GH.push_back( HyperBasevector( KX, p ) );
     vec<ref_loc> readlocs;
     long_logging_control log_control( ref, &readlocs, "", "" );
     long_heuristics heur( "" );
     heur.K2_FORCE = K2;
     String TMP = "tmp.xxx";
     LongHyper( "", correctede, cpartner, shb, heur, log_control,
          logc, TMP, USE_OLD_LRP_METHOD );
     PRINT( shb.EdgeObjectCount( ) );

     shb.DeleteReverseComplementComponents(logc);

     // Delete low coverage edges.

     if (DELETE_LOW_COV) 
     {    heur.LC_CAREFUL = True;
          shb.DeleteLowCoverage( heur, log_control, logc );    }

     // Trim hanging ends.

     const double junk_ratio = 10.0;
     const int max_del = 2500;
     shb.TrimHangingEnds( max_del, junk_ratio, heur, logc );
     
     // Clean up.

     shb.RemoveSmallComponents(5000);

     // Write assembly.

     PRINT( shb.EdgeObjectCount( ) );
     PRINT( shb.NComponents( ) );
     if ( shb.Acyclic( ) ) cout << "assembly is acyclic" << endl;
     Ofstream( out1, "assembly.dot" );
     shb.PrintSummaryDOT0w( out1, True, False, True );
     Ofstream( out2, "assembly.fasta" );
     for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
          shb.EdgeObject(e).Print( out2, e );
     BinaryWriter::writeFile( "assembly.shbv", shb );

     // Evaluate assembly.

     vec<int> inv( shb.EdgeObjectCount( ), -1 );
     int eval_verb = 0;
     int iAlignerK = 501;
     cout << Date( ) << ": begin evaluation" << endl;
     AssemblyError err = EvalAssembly( ref, shb, inv, cout, iAlignerK, eval_verb );

     // Summarize results.

     cout << "\nSUMMARY\n";
     ostringstream eout;
     eout << 1000 * double( shb.EdgeObjectCount( ) ) / genome[0].isize( );
     cout << "edges per kb = " << eout.str( ) << endl;
     cout << err.GetLeftGaps( ).size( ) + err.GetRightGaps( ).size( ) << " end gaps" 
          << endl;
     cout << err.GetGaps( ).size( ) << " internal gaps" << endl;
     int errs = err.GetIndels( ).isize( ) + err.GetSubs( );
     cout << "error rate = " << PERCENT_RATIO( 3, errs, genome[0].isize( ) ) << endl;
     cout << "total time used = " << TimeSince(clock) << endl;

     // Done.

     cout << "\n" << Date( ) << ": done\n" << endl;
     Scram(0);    }



// vim: sw=5
