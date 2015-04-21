///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Jul 31, 2014 - <crdhelp@broadinstitute.org>
//

#ifndef SCORES_H_
#define SCORES_H_
#include "paths/long/hic/CanonicalTestOrder.h"
#include "paths/long/hic/ScaffoldPile.h"

struct Scores {
    Scores() = default;
    Scores( size_t& ss, double ll, double lla, double chisq, size_t count = 0 ) :
        mSumScore(ss), mLogLikelihood(ll), mLogLikelihoodA(lla),
        mChiSq(chisq), mCount(count), mBestBit(0u) {};

    enum Metric { SUM, LOG_LIK, LOG_LIK_A, CHISQ };

    double operator()( Metric const m = SUM ) const {
        if ( m == SUM ) return mSumScore;
        else if ( m == LOG_LIK ) return mLogLikelihood;
        else if ( m == LOG_LIK_A ) return mLogLikelihoodA;
        else if ( m == CHISQ ) return mChiSq;
        else FatalErr("BUG: bad metric specified in Scores::operator()");
    }

    friend ostream& operator<<(ostream& out, Scores const& self) {
        out << "[sum " << ( ( self.mBestBit & 0x1 ) ? "*" : "" )
                << self.mSumScore << "]  [logLik " <<
                ( ( self.mBestBit & 0x2 ) ? "*" : "" )
                << self.mLogLikelihood <<
                "]  [1/x " <<
                ( ( self.mBestBit & 0x4 ) ? "*" : "" )
                << self.mLogLikelihoodA << "]" <<
                "  [chisq " << ((self.mBestBit & 0x8)?"*":"")
                << self.mChiSq << "] [count " << self.mCount << "]";
        return out;
    }

    size_t mSumScore;
    double mLogLikelihood;
    double mLogLikelihoodA;
    double mChiSq;
    size_t mCount;
    uchar mBestBit;
};

using VecScores = vec<Scores>;


Scores scoreHypothesis( vec<TaggedLineId> const& hyp,
        std::unordered_map<int,int>const& lineLengths,
        LineLocPairVec const& linePairs, SepHistogram const& model,
        int trace_id = -1 );

Scores testLinePairSpecific(
        std::pair<int,bool> configuration,
        TaggedLineId const& root, TaggedLineId const& test,
        unordered_map<int,int> const& lineLengths,
        LineLocPairVec const& linePairs, SepHistogram const& model);

VecScores testLinePair( TaggedLineId const& A, TaggedLineId const& B,
        std::unordered_map<int,int> const& lineLengths,
        LineLocPairVec const& linePairs, SepHistogram const& model);


pair<size_t,size_t> winningPairIndex( vec<Scores> scores );


pair<size_t,size_t> winningPairOrderOrient( vec<Scores> const scores,
        double& ratio, Scores::Metric metric = Scores::LOG_LIK_A );

pair<size_t,size_t> winningPairScore( vec<Scores> const scores,
        double& ratio, Scores::Metric metric = Scores::LOG_LIK_A );

double winningPairRatio( vec<Scores> scores, double lenSum, Scores::Metric metric );

double winningOrderRatio(vec<Scores> scores, Scores::Metric metric, Bool VERBOSE = False);


vec<TaggedLineId> scaffoldToHypothesis( vec<int> scaff, bool flip );

Scores testScaffPairSpecific( ScaffoldPile const& scaffPile,
        std::pair<int,bool> configuration, int root, int test,
        std::unordered_map<int,int> const& lineLengths,
        LineLocPairVec const& linePairs, SepHistogram const& model,
        LineInfoVec const& liv);

VecScores testScaffPair( ScaffoldPile const& scaffPile,
        int A, int B, unordered_map<int,int> const& lineLengths,
        LineLocPairVec const& linePairs, SepHistogram const& model,
        LineInfoVec const& liv);

#endif /* SCORES_H_ */
