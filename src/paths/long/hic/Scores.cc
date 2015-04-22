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

#include "paths/long/hic/HiCDefs.h"
#include "paths/long/hic/Scores.h"

double andersondarling(vec<size_t>& offsets, size_t hyp_len, int traceid = -1 )
{
    auto const badvalue = std::numeric_limits<double>::max();

    size_t low = 1000;
    size_t hi = hyp_len;

    if ( offsets.size() == 0 ) return badvalue;

    Sort(offsets);
    vec<long double> cumsum;
    vec<size_t> impulses;

    size_t i = 0;
    while (offsets[i] < low) i++;

    if ( i >= offsets.size() ) return badvalue;

    size_t last = offsets[i];
    long double sum = 0;
    for ( ; i < offsets.size(); ++i ) {
        auto off = offsets[i];
        if ( off != last ) {
            impulses.push_back(last);
            cumsum.push_back(sum);
        }
        sum += 1.;
        last = off;
    }
    impulses.push_back(offsets.back());
    cumsum.push_back(sum);

    for ( auto& c : cumsum )
        c /= sum;

    vec<long double> logx;
    for ( auto& i : impulses ) {
        if ( i <= low ) logx.push_back(log(1.01));
        else logx.push_back(log(i-low));
    }


    for ( auto& l : logx ) l /= logx.back();

    size_t n = cumsum.size() - 1;       // normalized so that the ends agree -- skip it
    long double ad = 0.;

    for ( size_t i = 0; i < n; ++i ) {
        ad += ( cumsum[i] - logx[i] )*( cumsum[i] - logx[i] )  /
                ( logx[i] * ( 1. - logx[i] ) ) ;
    }
    ad /= n;

    if ( traceid > -1 ) {
#pragma omp critical
        {
            Ofstream(out, "test-"+ToString(traceid)+".txt");
            for ( size_t i = 0; i < cumsum.size(); ++i )
                out << "CUMSUM " << impulses[i] << " " << cumsum[i] << endl;
            for ( size_t i = 0; i < logx.size(); ++i )
                out << "LOGX " << impulses[i] << " " << logx[i] << endl;
        }
    }

    return ad;
}

Scores scoreHypothesis( vec<TaggedLineId> const& hyp,
        std::unordered_map<int,int> const& lineLengths,
        LineLocPairVec const& linePairs, SepHistogram const& model,
        int trace_id )
{
    size_t hyp_len = std::accumulate(hyp.begin(), hyp.end(), 0u,
            [&lineLengths]( unsigned sum, TaggedLineId const& a )
            { return sum + lineLengths.at(a.getId()); } );
//    vec<size_t> raw_offsets;

    bool trace = false; //trace_id != -1;
    // build distribution given implied separations/adjacencies
//    SepHistogram seps;
//    seps.CloneShape(model);
    size_t count = 0;
    size_t sum_score = 0;
    double log_likelihood = 0.;
    double log_likelihood_a = 0.;
#pragma omp parallel for reduction(+:sum_score,count,log_likelihood,log_likelihood_a)
    for ( size_t i = 0; i < linePairs.size(); ++i ) {
        LineLocPair const& pair = linePairs[i];
        LineLoc ll1 = pair.ll1();
        LineLoc ll2 = pair.ll2();

        bool found = false;
        int offset = 0;
        if ( ll1.getLineId() == ll2.getLineId() ) {
            for ( auto itr1 = hyp.begin(); !found && itr1 != hyp.end(); ++itr1 ) {
                if ( itr1->getId() == ll1.getLineId() ) {
                    offset = std::abs( ll1.getLeftOffset() - ll2.getLeftOffset() );
//                    seps.Add( offset );
//                    if ( static_cast<unsigned>(offset) < hyp_len ) raw_offsets.push_back(offset);
                    if (trace) cout << "COUNT-" + ToString(trace_id) << " " << offset << endl;
                    count++;
                    found=true;
                }
            }
        } else {
            for ( auto itr1 = hyp.begin(); !found && itr1 != hyp.end(); ++itr1 ) {
                if ( itr1->getId() == ll2.getLineId() ) std::swap( ll1, ll2 );
                if ( itr1->getId() == ll1.getLineId() ) {
                    offset = itr1->isTag() ? ll1.getLeftOffset() : ll1.getRightOffset();
                    for ( auto itr2 = itr1+1; !found && itr2 != hyp.end(); ++itr2 ) {
                        if ( itr2->getId() == ll2.getLineId() ) {
                            // we found the pair -- add it
                            offset += itr2->isTag() ? ll2.getRightOffset() : ll2.getLeftOffset();
//                            seps.Add(offset);
//                            if ( static_cast<unsigned>(offset) < hyp_len ) raw_offsets.push_back(offset);
                            if (trace) cout << "COUNT-" + ToString(trace_id) << " " << offset << endl;
                            count++;
                            found = true;
                            break;
                        } else {
                            offset += lineLengths.at(itr2->getId());
                        }
                    }
                }
            }
        }

        if ( found && offset > 0 ) {
            sum_score += offset;
            double pval = model(offset);
            if ( pval > 1e-8 ) {
                log_likelihood += log( pval );
                if ( trace ) cout << "LOGLIK-"+ToString(trace_id) << " "
                    << offset << " " << -log(pval) << " "<< -log_likelihood << endl;
            }
            double pval_a = model.Analytic(offset);
            if ( pval_a > 1e-8 ) log_likelihood_a += log(pval_a );
        }
    }
//    seps.Smooth(4.);
//    if ( trace ) seps.Print(cout, "HIST-" +
//            ToString(trace_id) + " " );

//    seps.PrintStats(cout, "query");
//    model.PrintStats(cout, "model");

    // score distribution
    return Scores( sum_score, -log_likelihood, -log_likelihood_a, 0., count);
//    return Scores( sum_score, -log_likelihood, -log_likelihood_a,
//            andersondarling(raw_offsets, hyp_len, trace_id), count);
//            chisq(seps,model,6000,300000), count );
}




Scores testLinePairSpecific(
        std::pair<int,bool> configuration,
        TaggedLineId const& root, TaggedLineId const& test,
        unordered_map<int,int> const& lineLengths,
        LineLocPairVec const& linePairs, SepHistogram const& model)
{
    vec<TaggedInt> hyp;
    ForceAssert(configuration.first == -1 || configuration.first == 1 );
    if ( configuration.first == - 1) {
        hyp.push_back( test );
        if ( configuration.second ) hyp.back().flipTag();
        hyp.push_back( root );
    } else {
        hyp.push_back( root );
        hyp.push_back( test );
        if ( configuration.second ) hyp.back().flipTag();
    }

    return scoreHypothesis( hyp, lineLengths, linePairs, model );
}

VecScores testLinePair( TaggedLineId const& A, TaggedLineId const& B,
        std::unordered_map<int,int> const& lineLengths,
        LineLocPairVec const& linePairs, SepHistogram const& model)
{
    vec<Scores> scores;
    vec<std::pair<int,bool>> configs(4);

    // ORDER MUST MATCH pairIndexToConfig -- need to move this to a class
    for ( size_t i = 0; i < configs.size(); ++i )
        configs[i] = CanonicalTestOrder::pairIndexToPairConfig(i);

    for ( auto const& config : configs ) {
        scores.push_back( testLinePairSpecific( config,
                A, B, lineLengths, linePairs, model ) );
    }

    return scores;
}


pair<size_t,size_t> winningPairIndex( vec<Scores> scores )
{
    vec<int> idx(scores.size(), vec<int>::IDENTITY );
    SortSync( scores, idx, []( Scores const& s1, Scores const& s2 ) {
        return s1(Scores::LOG_LIK_A) < s2(Scores::LOG_LIK_A);
    });

    return make_pair(idx[0],idx[1]);
}


pair<size_t,size_t> winningPairOrderOrient( vec<Scores> const scores,
        double& ratio, Scores::Metric metric)
{
    auto idx = winningPairIndex(scores);

    int idx0 = idx.first;
    int idx1 = idx0 ^ 1u;       // 0<->1 and 2<->3

    auto s0 = scores[idx0](metric);
    auto s1 = scores[idx1](metric);

    if ( metric == Scores::LOG_LIK_A || metric == Scores::LOG_LIK ) {
        ratio = s0 - s1;
    } else {
        // NOTE**** we're not using this, but it's a different "sense"
        // (smaller is better) than some of the previous work.
//        ratio = scores[idx0](metric) / scores[idx1](metric);
        ratio = scores[idx0](metric) / scores[idx1](metric);
    }

    return make_pair(idx0,idx1);
}
#if 0
pair<size_t,size_t> winningPairOrderOrient( vec<Scores> const scores,
        double& ratio, Scores::Metric metric)
{
//    auto idx = winningPairIndex(scores);
//    PRINT2(scores[idx.first](metric),scores[idx.second](metric));
//    if ( scores[0].mCount > 0 &&
//            (scores[idx.second](metric) - scores[idx.first](metric)) / scores[idx.first](metric) < .002 ) {
//        ratio = 0.;
//        return idx;
//    }
    auto s0 = scores[0](metric);
    auto s1 = scores[1](metric);
    auto s2 = scores[2](metric);
    auto s3 = scores[3](metric);

    auto ab = s0 + s2;
    auto ba = s1 + s3;

    int idx0, idx1;
    if ( ab < ba ) {
        idx0 = idx1 = 0;
        if ( s0 < s2 ) idx1 = 2; else idx0 = 2;
    } else {
        idx0 = idx1 = 1;
        if ( s1 < s3 ) idx1 = 3; else idx0 = 3;
    }

    if ( metric == Scores::LOG_LIK_A || metric == Scores::LOG_LIK ) {
//        ratio = scores[idx0](metric) - scores[idx1](metric);
        if ( ab < ba ) ratio = ab - ba; else ratio = ba - ab;
    } else {
        // NOTE**** we're not using this, but it's a different "sense"
        // (smaller is better) than some of the previous work.
//        ratio = scores[idx0](metric) / scores[idx1](metric);
        ratio = scores[idx0](metric) / scores[idx1](metric);
    }

    return make_pair(idx0,idx1);
}
#endif

pair<size_t,size_t> winningPairScore( vec<Scores> const scores,
        double& ratio, Scores::Metric metric)
{
    auto idxpair = winningPairIndex( scores );
    size_t count = scores[0].mCount;

    if ( count == 0 ) ratio = 0.;
    else {

    double norm0 = 10.*scores[idxpair.first](metric) / (1.*count);
    double norm1 = 10.*scores[idxpair.second](metric) / (1.*count);

    ratio = norm0 - norm1;
    }

    return idxpair;
}


double winningPairRatio( vec<Scores> scores, double lenSum, Scores::Metric metric )
{
    auto idx = winningPairIndex(scores);

    if ( metric == Scores::LOG_LIK_A || metric == Scores::LOG_LIK ) {
        auto count = scores[idx.first].mCount;
        // scores are negative log, so this is equivalent to p(ALT) / p(WINNER)
        return exp(scores[idx.first](metric)-scores[idx.second](metric)) * count/100000;
    } else
            return static_cast<double>(scores[idx.second](metric))
                    / scores[idx.first](metric);
}

double winningOrderRatio(vec<Scores> scores, Scores::Metric metric, Bool VERBOSE)
{
    double ab = scores[0](metric) + scores[2](metric);
    double ba = scores[1](metric) + scores[3](metric);

    if ( VERBOSE ) {
        PRINT2(ab,ba);
    }

    auto idx = winningPairIndex(scores).first;

    if ( metric == Scores::LOG_LIK_A || metric == Scores::LOG_LIK ) {
       auto count =  scores[0].mCount;
       double val;
       if ( idx == 0 || idx == 2 ) val=ab - ba;
       else val=ba - ab;
       return val;
    } else {
        if ( idx == 0 || idx == 2 ) return (ba/ab); return (ab/ba);
    }
}


vec<TaggedLineId> scaffoldToHypothesis( vec<int> scaff, bool flip )
{
    vec<TaggedLineId> hyp;
    hyp.reserve(scaff.size());
    if ( flip ) {
        for ( auto itr = scaff.rbegin(); itr != scaff.rend(); ++itr )
            hyp.emplace_back( *itr, true );
    } else {
        for ( auto itr = scaff.begin(); itr != scaff.end(); ++itr )
            hyp.emplace_back( *itr, false );
    }

    return hyp;
}

vec<TaggedLineId> canonicalizeHypothesis(vec<TaggedLineId> hyp,
        LineInfoVec const& liv, bool trace = false )
{
    for ( auto& line : hyp ) {  // copy of hyp
        if ( liv[line.getId()].mIsRC ) {
            if ( trace ) cout << "canonicalizing " << line << " to ";
            line.setId(liv[line.getId()].mLineIdRC);
            if ( trace ) cout << line << " -> ";
            line.flipTag();
            if ( trace ) cout << line << endl;
        }
    }
    return hyp;
}


Scores testScaffPairSpecific(
        ScaffoldPile const& scaffPile,
        std::pair<int,bool> configuration,
        int root, int test,
        std::unordered_map<int,int> const& lineLengths,
        LineLocPairVec const& linePairs, SepHistogram const& model,
        LineInfoVec const& liv)
{
    vec<TaggedInt> hyp;
    ForceAssert(configuration.first == -1 || configuration.first == 1 );
    if ( configuration.first == -1) {   // test left of root
        hyp.append( scaffoldToHypothesis(scaffPile(test), configuration.second ) );
        hyp.append( scaffoldToHypothesis(scaffPile(root), false));
    } else {
        hyp.append( scaffoldToHypothesis(scaffPile(root), false ));
        hyp.append( scaffoldToHypothesis(scaffPile(test), configuration.second ));
    }

    return scoreHypothesis( canonicalizeHypothesis(hyp,liv), lineLengths, linePairs, model );
}

VecScores testScaffPair( ScaffoldPile const& scaffPile,
        int A, int B,
        unordered_map<int,int> const& lineLengths,
        LineLocPairVec const& linePairs, SepHistogram const& model,
        LineInfoVec const& liv)
{
    vec<Scores> scores;
    vec<std::pair<int,bool>> configs(4);

    // ORDER MUST MATCH pairIndexToConfig -- need to move this to a class
    for ( size_t i = 0; i < configs.size(); ++i )
        configs[i] = CanonicalTestOrder::pairIndexToPairConfig(i);

    for ( auto const& config : configs ) {
        scores.push_back( testScaffPairSpecific( scaffPile, config,
                A, B, lineLengths, linePairs, model, liv ) );
    }

    return scores;
}
