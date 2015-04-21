///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>
#include "MainTools.h"
#include "paths/long/hic/HiCDefs.h"
#include "VecUtilities.h"
#include <iosfwd>
#include "system/System.h"
#include <queue>
#include <unordered_set>
#include "paths/long/hic/TestTree.h"
#include "paths/long/hic/HiCExperiment.h"
#include "paths/long/hic/LineMapReader.h"

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
    Scores( size_t& ss, double ll, double lla ) :
        mSumScore(ss), mLogLikelihood(ll), mLogLikelihoodA(lla), mBestBit(0u) {};

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
    uchar mBestBit;
};

Scores scoreHypothesis( vec<TaggedLineId> const& hyp,
        std::unordered_map<int,int> lineLengths,
        LineLocPairVec const& linePairs, SepHistogram const& model,
        int trace_id = -1 )

{
    bool trace = trace_id != -1;
    size_t count = 0;
    size_t sum_score = 0;
    double log_likelihood = 0.;
    double log_likelihood_a = 0.;
    for ( LineLocPair const& pair : linePairs ) {
        LineLoc ll1 = pair.ll1();
        LineLoc ll2 = pair.ll2();

        int offset = 0;
        bool found = false;
        if ( ll1.getLineId() == ll2.getLineId() ) {
            offset = std::abs( ll1.getLeftOffset() - ll2.getLeftOffset() );
            if (trace) cout << "COUNT-" + ToString(trace_id) << " " << offset << endl;
            count++;
            found = true;
        } else {
            for ( auto itr1 = hyp.begin(); !found && itr1 != hyp.end(); ++itr1 ) {
                if ( itr1->getId() == ll2.getLineId() ) std::swap( ll1, ll2 );
                if ( itr1->getId() == ll1.getLineId() ) {
                    offset = itr1->isTag() ? ll1.getLeftOffset() : ll1.getRightOffset();
                    for ( auto itr2 = itr1+1; !found && itr2 != hyp.end(); ++itr2 ) {
                        if ( itr2->getId() == ll2.getLineId() ) {
                            // we found the pair -- add it
                            offset += itr2->isTag() ? ll2.getRightOffset() : ll2.getLeftOffset();
                            if (trace) cout << "COUNT-" + ToString(trace_id) << " " << offset << endl;
                            count++;
                            found = true;
                            break;
                        } else {
                            offset += lineLengths[itr2->getId()];
//                            if (trace) cout << "INTERMEDIATE!" << endl;
                        }
                    }
                }
            }
        }

//        if (!found && offset > 0 ) cout << "HOLY SHIT!" << endl;

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

    // score distribution
    return Scores( sum_score, -log_likelihood, -log_likelihood_a );
}

Scores testLinePairSpecific(
        std::pair<int,bool> configuration,
        TaggedLineId const& A, TaggedLineId const& B,
        std::unordered_map<int,int> lineLengths,
        LineLocPairVec const& linePairs, SepHistogram const& model, int trace_id = -1)
{
    vec<TaggedInt> hyp;
    ForceAssert(configuration.first == -1 || configuration.first == 1 );
    if ( configuration.first == - 1) {
        hyp.push_back( B );
        if ( configuration.second ) hyp.back().flipTag();
        hyp.push_back( A );
    } else {
        hyp.push_back( A );
        hyp.push_back( B );
        if ( configuration.second ) hyp.back().flipTag();
    }

    auto score = scoreHypothesis( hyp, lineLengths, linePairs, model, trace_id );
    return score;
}

//std::pair<int, bool> testLinePair( TaggedLineId const& A,
vec<Scores> testLinePair( TaggedLineId const& A,
        TaggedLineId const& B,
        std::unordered_map<int,int> lineLengths,
        LineLocPairVec const& linePairs, SepHistogram const& model)
{
    vec<Scores> scores;
    vec<std::pair<int,bool>> configs;

    // right forw
    configs.emplace_back( 1, false );
    // left forw
    configs.emplace_back( -1, false );
    // right rev
    configs.emplace_back( 1, true );
    // left rev
    configs.emplace_back( -1, true );

    int trace_id = -1;
    for ( auto const& config : configs ) {
        scores.push_back( testLinePairSpecific( config,
                A, B, lineLengths, linePairs, model, trace_id ) );
    }

    return scores;

#if 0
    vec<int> idx(scores.size(), vec<int>::IDENTITY );
    SortSync( scores, idx, []( Scores const& s1, Scores const& s2 ) {
        return s1.mLogLikelihood < s2.mLogLikelihood;
    });

    return configs[idx.front()];
#endif
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
    CommandArgument_IntSet_OrDefault(TRACE_LINE_IDS, "");
    CommandArgument_String_OrDefault_Doc(DUMPHEAD,"",
            "dump some debugging files with this head");
    CommandArgument_Double_OrDefault_Doc(SUBSAMPLE, 1.0, "subsample factor");
    EndCommandArguments;

    ForceAssertGe(SUBSAMPLE, 0.0);
    ForceAssertLe(SUBSAMPLE, 1.0);

    cout << Date() << ": reading a.s.lines.map" << endl;
    auto lineLoci = LineMapReader(ADIR + "/a.s.lines.map").LineLoci();

    HiCExperiment hic(
            ADIR + "/a.s.hbx",
            ADIR + "/a.s.inv",
            HICDIR + "/" + EdgeIdToLineIdMap::gFileName,
            HICDIR + "/" + LineInfo::gFileName,
            ADIR + "/" + "a.s.lines",
            HICDIR + "/hic.hitrates",
            SEPDIST
            );

    vec<LineMapReader::LineLocus> longLineLoci;
    for ( auto itr = lineLoci.begin(); itr != lineLoci.end()-1; ++itr ) {
        if ( itr[0].len > 100000 && itr[1].len > 100000 ) {
            longLineLoci.push_back(itr[0]);
            longLineLoci.push_back(itr[1]);
        }
    }

#if 0
    // make a list of lines > 100k
    std::copy_if( lineLoci.begin(), lineLoci.end(),
            std::back_inserter( longLineLoci ),
            [](LineMapReader::LineLocus const& l) {
                return ( l.len > 100000 ); } );
#endif

    cout << "there are " << longLineLoci.size() << " suprathreshold loci" << endl;

    vec<int> testLines;
    for ( auto const& lineLocus : longLineLoci )
        testLines.push_back(canonicalLine(lineLocus.lineId,hic.mLineInfoVec));

    PRINT(testLines.size());


    // pull hiCPairs incident on these lines
    cout << Date() << ": reading HiC pairs for selected lines" << endl;
    HICVec hiCPairs;
    getHiCPairs( HICPAIRS, hic.mEdgeToLineMap, &hiCPairs, true, testLines);
    cout << Date() << ": done reading " << hiCPairs.size() << " pairs." << endl;


    // turn pairs of the form:
    // ( (edge1, offset1), (edge2, offset2)
    // into:
    // ( (line1, loffset1, roffset1), (line2, loffset2, roffset2) )
    // including only contacts where both lines are in the list.
    cout << Date() << ": calculating edge pairs -> line pairs" << endl;
    LineOffsetCalc offCalc(hic.mHBX, hic.mLines, hic.mInv,
            hic.mEdgeToLineMap, hic.mRuler, LineOffsetCalc::THREAD_SAFE);

    cout << Date() << ": preloading cache of line lengths and edge<->line map" << endl;
    for ( auto const& testLine : testLines )
        offCalc.cachePreload(testLine);

    cout << Date() << ": scrub contacts not in our set of interest" << endl;
    vec<Bool> to_delete(hiCPairs.size(), false );
#pragma omp parallel for
    for ( size_t i = 0; i < hiCPairs.size(); ++i ) {
        int line1 = hic.mEdgeToLineMap[hiCPairs[i].el1().getEdgeId()];
        int line2 = hic.mEdgeToLineMap[hiCPairs[i].el2().getEdgeId()];
        if ( !Member(testLines, line1 ) || !Member(testLines, line2) )
            to_delete[i] = true;
    }
    EraseIf(hiCPairs, to_delete);

    cout << Date() << ": calculating..." << endl;
    if ( TRACE_LINE_IDS.size() && TRACE_LINE_IDS.size() != 2 ) FatalErr("should only trace 2 right now");

    std::ostringstream traceLog;
    LineLocPairVec linePairs(hiCPairs.size());
    std::map<int, triple<int,int,int> > traceCounts;
    size_t ab=0, ba=0, ab_=0, b_a=0;
#pragma omp parallel for
    for ( size_t i = 0; i < hiCPairs.size(); ++i ) {
        if ( omp_get_thread_num() == 0 && i % 100000 == 0) cout << i << " of " << hiCPairs.size() / omp_get_num_threads() << endl;
        auto const& hiCPair = hiCPairs[i];
        LineLocPair tmp( offCalc(hiCPair.el1()), offCalc(hiCPair.el2()) );
        linePairs[i]=tmp;

        // trace, if requested
        if ( TRACE_LINE_IDS.size() ) {
            auto lineId1 = tmp.ll1().getLineId();
            auto lineId2 = tmp.ll2().getLineId();
            if ( Member( TRACE_LINE_IDS, lineId1 ) &&
                    Member( TRACE_LINE_IDS, lineId2 ) ) {
#pragma omp critical
                {
                LineLoc A, B;
                if ( lineId1 == TRACE_LINE_IDS[0] ) {
                    ForceAssertEq(lineId2, TRACE_LINE_IDS[1]);
                    A = tmp.ll1();
                    B = tmp.ll2();
                } else  {
                    ForceAssertEq(lineId1, TRACE_LINE_IDS[1]);
                    ForceAssertEq(lineId2, TRACE_LINE_IDS[0]);
                    A = tmp.ll2();
                    B = tmp.ll1();
                }
                ab += A.getRightOffset() + B.getLeftOffset();
                ba += A.getLeftOffset() + B.getRightOffset();
                ab_ += A.getRightOffset() + B.getRightOffset();
                b_a += A.getLeftOffset() + B.getLeftOffset();
                triple<int,int,int> tmp1{0,0,0};
                if ( traceCounts.find( lineId1 ) != traceCounts.end() ) {
                    tmp1 = traceCounts[lineId1];
                }
                tmp1.first += tmp.ll1().getLeftOffset();
                tmp1.second += tmp.ll1().getRightOffset();
                tmp1.third += 1;
                traceCounts[lineId1] = tmp1;

                triple<int,int,int> tmp2{0,0,0};
                if ( traceCounts.find( lineId2 ) != traceCounts.end() ) {
                    tmp2 = traceCounts[lineId2];
                }
                tmp2.first += tmp.ll2().getLeftOffset();
                tmp2.second += tmp.ll2().getRightOffset();
                tmp2.third += 1;
                traceCounts[lineId2] = tmp2;
                }
            }
        }
    }
    cout << Date() << ": done" << endl;

    if ( TRACE_LINE_IDS.size() ) {
        for ( auto trace : traceCounts ) {
            auto lineId = trace.first;
            auto tmp = trace.second;
            cout << "TRACE: lineId=" << lineId << ", leftAvg=" <<
                    tmp.first / tmp.third << ", rightAvg=" <<
                    tmp.second / tmp.third << endl;
        }
        PRINT4(ab,ba,ab_,b_a);
    }

    cout << "linePairs.size()=" << linePairs.size() << endl;

    // length of the lines in question
    std::unordered_map<int, int> testLineLengths;
    for ( auto const lineId : testLines ) {
        testLineLengths[lineId] = GetLineLength(hic.mLines[lineId], hic.mRuler);
    }

    size_t numLoci = longLineLoci.size();
    // the following should really be numLoci * (numLoci-1) / 2, but then
    // addressing it would be prone to downstream errors.
    vec<vec<Scores>> all_results(numLoci * numLoci);
#pragma omp parallel for
    for (size_t i = 0; i < longLineLoci.size()-1; i+=2) {
        size_t j = i+1;
        auto const& locus1 = longLineLoci[i];
        auto const& locus2 = longLineLoci[j];

        if ( TRACE_LINE_IDS.size() &&
                ( !Member(TRACE_LINE_IDS,locus1.lineId) ||
                !Member(TRACE_LINE_IDS, locus2.lineId) ) ) continue;

        if ( locus1.chr != locus2.chr ||
                std::abs(locus1.start - locus2.start) > 2000000 )
            continue;

        if (SUBSAMPLE < 1.0 && randint(65536) > SUBSAMPLE*65536)
            continue;

        // for now, because these are taken from the map, they
        // should all be facing in one direction
        TaggedLineId line1(locus1.lineId, false);
        TaggedLineId line2(locus2.lineId, false);
        auto results = testLinePair(line1, line2, testLineLengths,
                linePairs, hic.mSepDist);
        all_results[i*numLoci+j] = results;
        if ( omp_get_thread_num() == 0 )
            cout << i << " of " << longLineLoci.size() / omp_get_num_threads()
            << endl;
    }

    auto idx2config = [](size_t idx)->String {
        if (idx==0) return "AB";
        else if ( idx == 1) return "BA";
        else if (idx == 2) return "AB'";
        else if ( idx==3) return "B'A";
        else FatalErr("weird index");
    };

    for (size_t i = 0; i < longLineLoci.size()-1; i+=2 ) {
        size_t j = i+1;
        auto const& locus1 = longLineLoci[i];
        auto const& locus2 = longLineLoci[j];
        auto results = all_results[i*numLoci+j];
        if ( results.size() == 0 ) continue;
        cout << "A=[" << locus1.lineId << "], B=[" << locus2.lineId << "]";
        for ( size_t i = 0; i < 4; ++i )
            cout << ", " << idx2config(i) << "=" <<
            results[i].mLogLikelihoodA;
        vec<int> idx(results.size(), vec<int>::IDENTITY );
        SortSync( results, idx, []( Scores const& s1, Scores const& s2 ) {
            return s1.mLogLikelihoodA < s2.mLogLikelihoodA;
        });
        cout << ", TOP=" << idx2config(idx.front());
        cout << ", ACTUAL=";
        if ( locus1.chr != locus2.chr ) cout << "n/a";
        else if ( locus1.start < locus2.start ) cout << "AB";
        else cout << "BA";
        cout << "  DIST=" << std::min(std::abs(locus1.start - locus2.stop),
                std::abs(locus1.stop - locus2.start));
        cout << endl;
    }

    return 0;
}
