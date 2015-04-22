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
#include <unordered_set>

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
    Scores( size_t& ss, double ll, double lla, double kld, size_t counts = 0 ) :
        mSumScore(ss), mLogLikelihood(ll), mLogLikelihoodA(lla),
        mKLDivergence(kld), mCounts(counts), mBestBit(0u) {};

    friend ostream& operator<<(ostream& out, Scores const& self) {
        out << "[sum " << ( ( self.mBestBit & 0x1 ) ? "*" : "" )
                << self.mSumScore << "]  [logLik " <<
                ( ( self.mBestBit & 0x2 ) ? "*" : "" )
                << self.mLogLikelihood <<
                "]  [1/x " <<
                ( ( self.mBestBit & 0x4 ) ? "*" : "" )
                << self.mLogLikelihoodA << "]";
        if ( self.mCounts > 0 ) {
                out << endl << "------> normalized by counts=" << self.mCounts <<
                "[sum " << ( ( self.mBestBit & 0x1 ) ? "*" : "" )
                << self.mSumScore/self.mCounts << "]  [logLik " <<
                ( ( self.mBestBit & 0x2 ) ? "*" : "" )
                << self.mLogLikelihood/self.mCounts <<
                "]  [1/x " <<
                ( ( self.mBestBit & 0x4 ) ? "*" : "" )
                << self.mLogLikelihoodA/self.mCounts << "]";
        } else  {
            cout << "NO EVIDENCE" << endl;
        }
        return out;
    }

    size_t mSumScore;
    double mLogLikelihood;
    double mLogLikelihoodA;
    double mKLDivergence;
    size_t mCounts;
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

        int offset = 0;
        bool found = false;
        if ( ll1.getLineId() == ll2.getLineId() ) {
            offset = std::abs( ll1.getLeftOffset() - ll2.getLeftOffset() );
            seps.Add( offset );
            if (trace) cout << "COUNT-" + ToString(trace_id) << " " << offset << endl;
            count++;
            found = true;
//            if ( offset > 400000 ) { cout << "WARNING!: "; PRINT3(ll1, ll2, offset); }
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
//                            if ( offset <= 10000 ) count_below_10k++; else count_above_10k++;
                            found = true;
                            break;
                        } else {
                            offset += lineLengths[itr2->getId()];
                        }
                    }
                }
            }
//            if ( done && offset > 400000 ) { cout << "WARNING!!!: "; PRINT3(ll1, ll2, offset); }
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
        } else if ( found && offset == 0 ) {
            if ( ll1.getLineId() != ll2.getLineId() ) cout << "NOVEL ";
            cout << "ZERO: " << ll1.getLineId() << "," << ll2.getLineId() << endl;
        }
    }
//    PRINT3(count,count_below_10k,count_above_10k);
    seps.Smooth(4.);
    if ( trace ) seps.Print(cout, "HIST-" +
            ToString(trace_id) + " " );

//    seps.PrintStats(cout, "query");
//    model.PrintStats(cout, "model");

    // score distribution
    return Scores( sum_score, -log_likelihood, -log_likelihood_a,
            relative_entropy(seps,model,0,1000000,trace_id), count );
}

};

int main( int argc, char* argv[] )
{
    RunTime( );

    BeginCommandArguments;
    CommandArgument_String_Doc(ADIR,"assembly dir.");
    CommandArgument_String_Doc(HICDIR,"hic.* files from CalcHiCLineHits");
    CommandArgument_IntSet_OrDefault_Doc(LINES,"{}",
            "line ids to test as potential neighbors");
    CommandArgument_StringSet_OrDefault_Doc(HYPOTHESES,"{}",
            "e.g. {ABCD,BC'DA'}");
    CommandArgument_IntSet_OrDefault_Doc(TRACE,"{}",
            "Hypotheses to trace e.g. {0,1,2}");
    CommandArgument_String_Doc(HICPAIRS,"Hi-C edge-pairing file.");
    CommandArgument_String_Doc(SEPDIST,"hic.*.sepdist file");
    CommandArgument_String_OrDefault_Doc(DUMPHEAD,"",
            "dump some debugging files with this head");
    CommandArgument_Bool_OrDefault_Doc(FILTERSELF,False,
            "filter out self-line hits");
    EndCommandArguments;

    cout << Date() << ": reading involution" << endl;
    vec<int> inv;
    BinaryReader::readFile(ADIR+"/a.inv", &inv);

    cout << Date() << ": reading edge to line id map" << endl;
    EdgeIdToLineIdMap edgeToLineMap;
    BinaryReader::readFile( HICDIR+"/"+EdgeIdToLineIdMap::gFileName,
            &edgeToLineMap );

    cout << Date() << ": reading lineInfo" << endl;
    LineInfoVec lineInfoVec;
    BinaryReader::readFile( HICDIR+"/"+LineInfo::gFileName,
            &lineInfoVec );

    cout << Date() << ": reading HyperBasevectorX(treme)" << endl;
    HyperBasevectorX hbx;
    BinaryReader::readFile( ADIR + "/a.hbx", &hbx );
    HBVXKmerRuler ruler(hbx);

    cout << Date() << ": reading lines"  << endl;
    LineVec lines;
    BinaryReader::readFile(ADIR + "/a.lines" ,&lines);

    cout << Date() << ": reading distribution" << endl;
    SepHistogram sepdist;
    BinaryReader::readFile(SEPDIST, &sepdist);
    sepdist.Smooth(4.);
    sepdist.PrintStats(cout);
    sepdist.FixAnalyticSum();
    if ( TRACE.size() ) sepdist.Print(cout, "HIST-M ");

    vec<int> testLines;
    for ( auto const line : LINES )
        testLines.push_back(canonicalLine(line,lineInfoVec));

    // pull hiCPairs incident on these lines
    cout << Date() << ": reading HiC pairs for selected lines" << endl;
    HICVec hiCPairs;
    getHiCPairs( HICPAIRS, edgeToLineMap, &hiCPairs, false, testLines);

    // turn pairs of the form:
    // ( (edge1, offset1), (edge2, offset2)
    // into:
    // ( (line1, loffset1, roffset1), (line2, loffset2, roffset2) )
    // including only contacts where both lines are in the list.
    LineOffsetCalc offCalc(hbx, lines, inv, edgeToLineMap, ruler);
    LineLocPairVec linePairs;
    for ( auto const& hiCPair : hiCPairs ) {
        LineLocPair tmp( offCalc(hiCPair.el1()), offCalc(hiCPair.el2()) );
        if ( !Member(testLines, tmp.first.getLineId()) ||
                !Member(testLines, tmp.second.getLineId()) ) continue;
        if ( tmp.first.getLineId() == tmp.second.getLineId() && FILTERSELF )
            continue;
        linePairs.push_back(tmp);
    }

    if ( DUMPHEAD != "" )
        BinaryWriter::writeFile(DUMPHEAD+".llpairs", linePairs);

    cout << "linePairs.size()=" << linePairs.size() << endl;

    // length of the lines in question
    cout << "Line Lengths: " << endl;
    std::unordered_map<int, int> testLineLengths;
    for ( auto const lineId : testLines ) {
        testLineLengths[lineId] = GetLineLength(lines[lineId], ruler);
        cout << "lineId " << lineId << ": " << testLineLengths[lineId] << endl;
    }

    if ( HYPOTHESES.size() > 0 ) {
        // score each hypothesis
        vec<Scores> hscores;
        for ( int i = 0; i < HYPOTHESES.isize(); ++i ) {
            vec<TaggedLineId> hyp = getHypothesis( HYPOTHESES[i], testLines );
            int trace = Member(TRACE,i) ? i : -1;
            auto score = scoreHypothesis( hyp, testLineLengths, linePairs, sepdist, trace );
            hscores.push_back(score);
        }

        for ( int i = 0; i < HYPOTHESES.isize(); ++i ) {
            vec<TaggedLineId> hyp = getHypothesis( HYPOTHESES[i], testLines );
            cout << "HYPOTHESIS " << i << " = " << HYPOTHESES[i] <<
                    " -> ";
            for ( auto const& id : hyp )
                cout << id << " ";
            cout << endl << " --> score: " << hscores[i] ;
            cout << endl;
        }
    } else {
        // no hypotheses specified
        // we assume ABCDEF... and score each pairwise transposition
        cout << "scoring possible transpositions: " << endl;
        String result_sum, result_lik, result_lika, result_kld;
        vec<TaggedLineId> hyp;
        for ( size_t inner = 0; inner < testLines.size(); ++inner ) {
            hyp.emplace_back( testLines[inner], false );
        }
//        for ( auto h : hyp ) cout << h.getId() << " ";
//        cout << endl;
        auto best_score = scoreHypothesis( hyp, testLineLengths, linePairs, sepdist );

        // scoring each pairwise transposition
        for ( size_t first = 0; first < testLines.size()-1; ++first ) {
            hyp.clear();
            for ( size_t inner = 0; inner < testLines.size(); ++inner ) {
                if ( inner == first ) hyp.emplace_back( testLines[inner+1], false );
                hyp.emplace_back( testLines[inner], false);
                if ( inner == first ) ++inner;
            }

//            for ( auto h : hyp ) cout << h.getId() << " ";
//            cout << endl;
            auto score = scoreHypothesis( hyp, testLineLengths, linePairs, sepdist );
            if ( score.mSumScore < best_score.mSumScore ) result_sum+="-"; else result_sum+="+";
            if ( score.mLogLikelihood < best_score.mLogLikelihood ) result_lik += "-"; else result_lik += "+";
            if ( score.mLogLikelihoodA < best_score.mLogLikelihoodA ) result_lika += "-"; else result_lika += "+";
            if ( score.mKLDivergence < best_score.mKLDivergence ) result_kld += "-"; else result_kld += "+";
        }

        cout << "trans based on SumScore: " << result_sum << endl;
        cout << "trans based on LogLikel: " << result_lik << endl;
        cout << "trans based on 1/x     : " << result_lika << endl;
//        cout << "based on KLDiverg: " << result_kld << endl;
        cout << endl;

        cout << "testing by flipping each line:" << endl;
        // now we score each single flip
        result_sum=result_lik=result_lika=result_kld="";
        for ( size_t i = 0; i < hyp.size(); ++i ) {
            vec<TaggedLineId> tmp(hyp);
            tmp[i].setTag(true);
            auto score = scoreHypothesis( tmp, testLineLengths, linePairs, sepdist );
            if ( score.mSumScore < best_score.mSumScore ) result_sum+="-"; else result_sum+="+";
            if ( score.mLogLikelihood < best_score.mLogLikelihood ) result_lik += "-"; else result_lik += "+";
            if ( score.mLogLikelihoodA < best_score.mLogLikelihoodA ) result_lika += "-"; else result_lika += "+";
            if ( score.mKLDivergence < best_score.mKLDivergence ) result_kld += "-"; else result_kld += "+";
        }
        cout << "flip based on SumScore: " << result_sum << endl;
        cout << "flip based on LogLikel: " << result_lik << endl;
        cout << "flip based on 1/x     : " << result_lika << endl;
//        cout << "based on KLDiverg: " << result_kld << endl;
        cout << endl;


    }

    cout << Date() << ": done" << endl;
}
