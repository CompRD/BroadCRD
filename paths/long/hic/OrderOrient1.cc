///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "paths/long/hic/HiCDefs.h"
#include "VecUtilities.h"
#include <iosfwd>
#include "system/System.h"
#include <queue>
#include <unordered_set>
#include "paths/long/hic/TestTree.h"
#include "paths/long/hic/HiCExperiment.h"

namespace {

std::ostream& print_orientation( std::ostream& out, int id1, int id2, int variant )
{
        String line1id("---" + ToString(id1) + "---");
        String line2id("---" + ToString(id2) + "---");
        if ( variant == 0 ) out << "<" + line1id << line2id + ">";
        else if (variant == 1 ) out << "<" + line1id << "<" + line2id;
        else if (variant == 2 ) out << line1id + ">" << line2id + ">";
        else if ( variant == 3 ) out << line1id + ">" << "<" + line2id;
        else FatalErr("bad variant passed");
        return out;
}


vec<TaggedLineId> getHypothesis( String symbolic, vec<int> testLines )
{
    vec<TaggedLineId> ret;

    for ( auto itr = symbolic.begin(); itr != symbolic.end();  ) {
        int idx = *itr - 'A';
        if ( idx >= testLines.isize() || idx < 0 ) {
            std::ostringstream s;
            s << "hypothesis " << symbolic << " contains symbols outside " <<
                    "of A.." << static_cast<char>(testLines.size()-1+'A') ;
            FatalErr(s.str());
        }
        ++itr;
        if ( itr != symbolic.end() && *itr == '\'' ) {
            ret.push_back( TaggedLineId(testLines[idx],true)  );
            ++itr;
        } else
            ret.push_back(TaggedLineId(testLines[idx], false) );
    }

    return ret;
}


struct Scores {
    Scores() = default;
    Scores( size_t& ss, double ll, double lla, double kld ) :
        mSumScore(ss), mLogLikelihood(ll), mLogLikelihoodA(lla), mKLDivergence(kld), mBestBit(0u) {};

    friend ostream& operator<<(ostream& out, Scores const& self) {
        out << "[sum " << ( ( self.mBestBit & 0x1 ) ? "*" : "" )
                << self.mSumScore << "]  [logLik " <<
                ( ( self.mBestBit & 0x2 ) ? "*" : "" )
                << self.mLogLikelihood <<
                "]  [1/x " <<
                ( ( self.mBestBit & 0x4 ) ? "*" : "" )
                << self.mLogLikelihoodA << "]";
        return out;
    }

    size_t mSumScore;
    double mLogLikelihood;
    double mLogLikelihoodA;
    double mKLDivergence;
    uchar mBestBit;
};

Scores scoreHypothesis( vec<TaggedLineId> const& hyp,
        std::unordered_map<int,int> lineLengths,
        LineLocPairVec const& linePairs, SepHistogram const& model,
        int trace_id = -1 )
{
    bool trace = trace_id != -1;
    // build distribution given implied separations/adjacencies
    SepHistogram seps;
    seps.CloneShape(model);
    size_t count = 0;
    size_t sum_score = 0;
    double log_likelihood = 0.;
    double log_likelihood_a = 0.;
    for ( LineLocPair const& pair : linePairs ) {
        LineLoc ll1 = pair.ll1();
        LineLoc ll2 = pair.ll2();

        bool found = false;
        int offset = 0;
        if ( ll1.getLineId() == ll2.getLineId() ) {
            offset = std::abs( ll1.getLeftOffset() - ll2.getLeftOffset() );
            seps.Add( offset );
            if (trace) cout << "COUNT-" + ToString(trace_id) << " " << offset << endl;
            count++;
            found=true;
        } else {
            for ( auto itr1 = hyp.begin(); !found && itr1 != hyp.end(); ++itr1 ) {
                if ( itr1->getId() == ll2.getLineId() ) std::swap( ll1, ll2 );
                if ( itr1->getId() == ll1.getLineId() ) {
                    offset = itr1->isTag() ? ll1.getLeftOffset() : ll1.getRightOffset();
                    for ( auto itr2 = itr1+1; !found && itr2 != hyp.end(); ++itr2 ) {
                        if ( itr2->getId() == ll2.getLineId() ) {
                            // we found the pair -- add it
                            offset += itr2->isTag() ? ll2.getRightOffset() : ll2.getLeftOffset();
                            seps.Add(offset);
                            if (trace) cout << "COUNT-" + ToString(trace_id) << " " << offset << endl;
                            count++;
                            found = true;
                            break;
                        } else {
                            offset += lineLengths[itr2->getId()];
                        }
                    }
                }
            }
        }

        if ( found ) {
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
    seps.Smooth(4.);
    if ( trace ) seps.Print(cout, "HIST-" +
            ToString(trace_id) + " " );

//    seps.PrintStats(cout, "query");
//    model.PrintStats(cout, "model");

    // score distribution
    return Scores( sum_score, -log_likelihood, -log_likelihood_a,
            relative_entropy(seps,model,0,1000000,trace_id) );
}

Scores testLinePairSpecific(
        std::pair<int,bool> configuration,
        TaggedLineId const& root, TaggedLineId const& test,
        std::unordered_map<int,int> lineLengths,
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

std::pair<int, bool> testLinePair( TaggedLineId const& root,
        TaggedLineId const& test,
        std::unordered_map<int,int> lineLengths,
        LineLocPairVec const& linePairs, SepHistogram const& model)
{
    vec<Scores> scores;
    vec<std::pair<int,bool>> configs;

    // left forw
    configs.emplace_back( -1, false );
    // left rev
    configs.emplace_back( -1, true );
    // right forw
    configs.emplace_back( 1, false );
    // right rev
    configs.emplace_back( 1, true );

    for ( auto const& config : configs ) {
        scores.push_back( testLinePairSpecific( config,
                root, test, lineLengths, linePairs, model ) );
    }

    vec<int> idx(scores.size(), vec<int>::IDENTITY );
    SortSync( scores, idx, []( Scores const& s1, Scores const& s2 ) {
        return s1.mLogLikelihood < s2.mLogLikelihood;
    });

    return configs[idx.front()];
}



};

int main( int argc, char* argv[] )
{
    RunTime( );

    BeginCommandArguments;
    CommandArgument_String_Doc(ADIR,"assembly dir.");
    CommandArgument_String_Doc(HICDIR,"hic.* files from CalcHiCLineHits");
    CommandArgument_String_Doc(HICPAIRS,"Hi-C edge-pairing file.");
    CommandArgument_String_Doc(SEPDIST,"hic.*.sepdist file");
    CommandArgument_Int_Doc(SEED,"seed line");
    CommandArgument_String_OrDefault_Doc(DUMPHEAD,"",
            "dump some debugging files with this head");
    EndCommandArguments;

    HiCExperiment hic(
            ADIR + "/a.s.hbx",
            ADIR + "/a.s.inv",
            HICDIR + "/" + EdgeIdToLineIdMap::gFileName,
            HICDIR + "/" + LineInfo::gFileName,
            ADIR + "/" + "a.s.lines",
            HICDIR + "/hic.hitrates",
            SEPDIST
            );

    cout << Date() << ": calculating cohort based on seed" << endl;
    auto cseed = canonicalLine(SEED, hic.mLineInfoVec);
    if ( cseed != SEED ) {
        cout << "INFO: seed canonicalized to " << cseed << endl;
    }
    vec<int> testLines;
    cout << "HitRateVec[" << cseed << "]=";
    auto rates = hic.mHitRates[cseed];
    if ( rates.size() == 0 ) {
        cout << "nothing to do!" << endl;
        return 0;
    }
    auto rateThresh = rates[0].mHitsPerKb * 0.10;
    for ( auto const& hit : rates ) {
        if ( hit.mHitsPerKb >= rateThresh ) {
                cout << hit << ", ";
                testLines.push_back( hit.mLineId );
        }
    }
    cout << endl;

    // checking whether cohort is canonical
    for ( auto const& test : testLines )
        cout << test << " => " << canonicalLine(test,hic.mLineInfoVec) <<
        "|~" << hic.mLineInfoVec[test].mLineIdRC << endl;

    // pull hiCPairs incident on these lines
    cout << Date() << ": reading HiC pairs for selected lines" << endl;
    HICVec hiCPairs;
    getHiCPairs( HICPAIRS, hic.mEdgeToLineMap, &hiCPairs, true, testLines);

    cout << Date() << ": done" << endl;
    LineOffsetCalc offCalc( hic.mHBX, hic.mLines, hic.mInv,
            hic.mEdgeToLineMap, hic.mRuler );
    LineLocPairVec linePairs;
    for ( auto const& hiCPair : hiCPairs ) {
        LineLocPair tmp( offCalc(hiCPair.el1()), offCalc(hiCPair.el2()) );
        if ( !Member(testLines, tmp.first.getLineId()) ||
                !Member(testLines, tmp.second.getLineId()) ) continue;
        linePairs.push_back(tmp);
    }

    if ( DUMPHEAD != "" )
        BinaryWriter::writeFile(DUMPHEAD+".llpairs", linePairs);

    cout << "linePairs.size()=" << linePairs.size() << endl;

    // length of the lines in question
    cout << "Line Lengths: " << endl;
    std::unordered_map<int, int> testLineLengths;
    for ( auto const lineId : testLines ) {
        testLineLengths[lineId] = GetLineLength(hic.mLines[lineId], hic.mRuler);
          cout << "lineId " << lineId << ": " << testLineLengths[lineId] << endl;
      }

    vec<int> stack( testLines.rbegin(), testLines.rend() );

    // rooted at the seed, make a tree of order and orientations
    Tree<TaggedLineId> root(TaggedLineId(cseed, false));

    PRINT2(cseed, root.Node() );

    while ( stack.size() ) {
        // pop the next top-scoring value off of the stack
        int testLine = stack.back();
        stack.pop_back();

        // use a binary search to figure out where to insert
        auto const& sepdist = hic.mSepDist;
        root.BinaryWalk([testLine,&testLineLengths,&linePairs,&sepdist](Tree<TaggedLineId>* t) -> int {
            auto res = testLinePair( t->Node(), TaggedLineId(testLine, false), testLineLengths, linePairs, sepdist);
            ForceAssert( res.first == -1 || res.first == 1 );
            if ( res.first == -1 ) {
                if ( t->isLeft() ) return res.first;
                else t->AddLeft( TaggedLineId( testLine, res.second ) );
                cout << "adding " << testLine << " left of " << t->Node() << endl;
            } else if ( res.first == 1 ) {
                if ( t->isRight() ) return res.first;
                else t->AddRight( TaggedLineId( testLine, res.second ));
                cout << "adding " << testLine << " right of " << t->Node() << endl;
            }
            return 0;
        } );
    }

    cout << "traversal:" << endl;
    root.DFWalk([](Tree<TaggedLineId>* t)->bool{
        cout << t->Node() << endl;
        return false;
    }, Tree<TaggedLineId>::Order::MIDDLE);


    return 0;
}
