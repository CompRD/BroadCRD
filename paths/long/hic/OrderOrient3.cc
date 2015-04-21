///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
//
// MakeDepend: archived

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
#include "paths/long/hic/Scores.h"
#include "graph/Digraph.h"




struct EdgeEvidence {
    double ratio;
    size_t count;
    size_t alt_config;
};

extern template class digraphVE<TaggedLineId,EdgeEvidence>;
using GraphType = digraphVE<TaggedLineId,EdgeEvidence>;

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

VecScores testLinePair( TaggedLineId const& A, TaggedLineId const& B,
        std::unordered_map<int,int> lineLengths,
        LineLocPairVec const& linePairs, SepHistogram const& model)
{
    vec<Scores> scores;
    vec<std::pair<int,bool>> configs;

    // ORDER MUST MATCH pairIndexToConfig -- need to move this to a class
    // AB
    configs.emplace_back( 1, false );
    // BA
    configs.emplace_back( -1, false );
    // AB'
    configs.emplace_back( 1, true );
    // B'A
    configs.emplace_back( -1, true );

    for ( auto const& config : configs ) {
        scores.push_back( testLinePairSpecific( config,
                A, B, lineLengths, linePairs, model ) );
    }

    return scores;
}

// order below MUST MATCH testLinePair -- this
// needs to be a class
String pairIndexToConfigString( size_t idx )
{
    if (idx==0) return "AB";
    else if ( idx == 1) return "BA";
    else if (idx == 2) return "AB'";
    else if ( idx==3) return "B'A";
    else FatalErr("weird index");
}

pair<TaggedLineId,TaggedLineId> pairIndexToConfig( size_t idx, TaggedLineId const& A, TaggedLineId const& B)
{
    auto Bpr = B;
    Bpr.flipTag();
    if ( idx == 0 ) return make_pair(A,B);
    else if ( idx == 1 ) return make_pair(B,A);
    else if ( idx == 2 ) return make_pair(A,Bpr);
    else if ( idx == 3 ) return make_pair(Bpr,A);
    else FatalErr("weird index");
}

pair<size_t,size_t> winningPairIndex( vec<Scores> scores )
{
    vec<int> idx(scores.size(), vec<int>::IDENTITY );
    SortSync( scores, idx, []( Scores const& s1, Scores const& s2 ) {
        return s1(Scores::LOG_LIK_A) < s2(Scores::LOG_LIK_A);
    });

    return make_pair(idx[0],idx[1]);
}

double winningPairRatio( vec<Scores> scores, Scores::Metric metric )
{
    auto idx = winningPairIndex(scores);

    if ( metric == Scores::LOG_LIK_A || metric == Scores::LOG_LIK ) {
        auto count = scores[idx.first].mCount;
        return exp( scores[idx.first](metric) -scores[idx.second](metric) )*count/100000.;
    } else
            return static_cast<double>(scores[idx.first](metric))
                    / scores[idx.second](metric);
}


//void RemoveTransitiveLinksFrom( digraphV<TaggedLineId>& graph, int v ) {
//    graph.From(v)
//}

void dumpGraph( GraphType const& graph, String const& filename, int cseed = -1, Bool singletons = False )
{
    Ofstream(OUT, filename);

    vec<String> vertex_labels(graph.N());
    vec<vec<String>> edge_labels(graph.N());
    vec<String> vertex_colors;
    for ( int v = 0; v < graph.N(); ++v ) {
        if ( graph.FromSize(v) != 0 || graph.ToSize(v) != 0 ) {
                vertex_labels[v] = ToString(graph.Vert(v).getId());
                if ( graph.Vert(v).isTag() ) vertex_labels[v] += "'";
        }
        for ( int j = 0; j < graph.FromSize(v); ++j ) {
            auto obj = graph.EdgeObjectByIndexFrom(v,j);
            edge_labels[v].push_back( ToString(obj.count) + " (" +
                    ToString(obj.ratio) + "|" +
                    pairIndexToConfigString(obj.alt_config) + ")" );
        }
    }

    if ( cseed >= 0 ) {
        for ( int v = 0; v < graph.N(); ++v ) {
            auto taggedLine = graph.Vert(v);
            vertex_colors.push_back( taggedLine.getId() == cseed &&
                    !taggedLine.isTag() ? "green" : "black" );
        }
    }
    graph.DOT_vl(OUT,vertex_labels,"", vec<vec<String>>(),
            vec<String>(), edge_labels, vertex_colors);
}

std::pair<bool,bool> inferLeftRight( int cseed, GraphType const& graph,
        TaggedLineId& left, TaggedLineId& right )
{
    double thresh=1e-5;

    ForceAssertGt(graph.N(),0);
    auto seed = graph.Vert(0);
    bool found_left = false, found_right = false;

    // list left nodes and confidences
    vec<int> lefty_idx;
    vec<double> lefty_vals;
    for ( int j = 0; j < graph.ToSize(0); ++j ) {
        auto ratio = graph.EdgeObjectByIndexTo(0,j).ratio;
        if ( ratio < thresh ) {
            lefty_vals.push_back( ratio );
            lefty_idx.push_back(graph.To(0)[j]);
        }
    }
    SortSync(lefty_vals, lefty_idx );

    // find winner on the left
    int win = 0;
    for ( ; win < lefty_idx.isize(); ++win ) {
        if ( !Meet2(graph.From(win),lefty_idx) ) break;
    }
    if ( win < lefty_idx.isize() ) {
        found_left = true;
        left = graph.Vert(lefty_idx[win]);
    }

    // list right nodes and confidences
    vec<int> righty_idx;
    vec<double> righty_vals;
    for ( int j = 0; j < graph.FromSize(0); ++j ) {
        auto ratio = graph.EdgeObjectByIndexFrom(0,j).ratio;
        if ( ratio < thresh ) {
            righty_vals.push_back( ratio );
            righty_idx.push_back(graph.From(0)[j]);
        }
    }
    SortSync(righty_vals, righty_idx );

    // find winner on the right
    win = 0;
    for ( ; win < righty_idx.isize(); ++win ) {
        if ( !Meet2(graph.To(win),righty_idx) ) break;
    }
    if ( win < righty_idx.isize() ) {
        found_right = true;
        right = graph.Vert(righty_idx[win]);
    }



    return make_pair(found_left, found_right );
}

size_t pathsBetween( GraphType const& graph, int v, int w )
{
    std::unordered_set<int> seen;
    std::vector<int> stack;
    stack.push_back(v);
    size_t hits = 0;

    while (stack.size()) {
        v = stack.back();
        stack.pop_back();
        if ( v == w ) {
            hits++;
        } else if ( seen.find(v) == seen.end() ) {
            seen.insert(v);
            for ( auto const from: graph.From(v) )
                stack.push_back(from);
        }
    }

    return hits;
}

std::pair<bool,bool> inferLeftRight2( int cseed, GraphType const& graph0,
        TaggedLineId& left, TaggedLineId& right )
{
    double thresh=1e-10;


    ForceAssertGt(graph0.N(),0);
    auto seed = graph0.Vert(0);
    bool found_left = false, found_right = false;

    GraphType graph(graph0);
    for (int v = 0; v < graph.N(); ++v ) {
        vec<int> to_delete;
        for (int j = 0; j < graph.FromSize(v); ++j )
            if ( graph.EdgeObjectByIndexFrom(v,j).ratio >= thresh )
                to_delete.push_back(j);
        graph.DeleteEdgesFrom(v, to_delete);
    }
    dumpGraph(graph, "graphTest.dot", cseed);


    // list left nodes and confidences
    vec<int> lefty_idx;
    vec<double> lefty_vals;
    vec<int> lefty_vert;
    for ( int j = 0; j < graph.ToSize(0); ++j ) {
        auto ratio = graph.EdgeObjectByIndexTo(0,j).ratio;
        if ( ratio < thresh ) {
            lefty_vals.push_back( ratio );
            lefty_idx.push_back(j);
        }
        lefty_vert.push_back(graph.To(0)[j]);
    }

    // find winner on the left
    if ( lefty_idx.size() > 0 ) {
        vec<int> lefty_betweens;
        for ( size_t i = 0; i < lefty_idx.size(); ++i ) {
            lefty_betweens.push_back( pathsBetween(graph, lefty_vert[lefty_idx[i]], 0) );
        }
        SortSync(lefty_betweens, lefty_idx, lefty_vals);
        for ( size_t i = 0; i < lefty_idx.size(); ++i ) {
            cout << "lefty " << i << " betweens=" << lefty_betweens[i] <<
                    ", vert=" << graph.Vert(lefty_vert[lefty_idx[i]]) << ", vals=" <<
                    lefty_vals[i] << endl;
        }
        // break ties
        int j = 1;
        while ( j < lefty_betweens.isize() && lefty_betweens[j] == lefty_betweens[0] ) j++;
        if ( j < lefty_betweens.isize() ) {
            lefty_vals.erase( lefty_vals.begin()+j, lefty_vals.end() );
            lefty_idx.erase( lefty_idx.begin()+j, lefty_idx.end() );
        }
        SortSync(lefty_vals, lefty_idx );
        for ( size_t i = 0; i < lefty_idx.size(); ++i ) {
            cout << "lefty " << i << " betweens=" << lefty_betweens[i] <<
                    ", vert=" << graph.Vert(lefty_vert[lefty_idx[i]]) << ", vals=" <<
                    lefty_vals[i] << endl;
        }

        left = graph.Vert(lefty_vert[lefty_idx[0]]);
        found_left = true;
    }


    // list right nodes and confidences
    vec<int> righty_idx;
    vec<double> righty_vals;
    vec<int> righty_vert;
    for ( int j = 0; j < graph.FromSize(0); ++j ) {
        auto ratio = graph.EdgeObjectByIndexFrom(0,j).ratio;
        if ( ratio < thresh ) {
            righty_vals.push_back( ratio );
            righty_idx.push_back(j);
        }
        righty_vert.push_back(graph.From(0)[j]);
    }

    // find winner on the right
    if ( righty_idx.size() > 0 ) {
        vec<int> righty_betweens;
        for ( size_t i = 0; i < righty_idx.size(); ++i ) {
            righty_betweens.push_back( pathsBetween(graph,0,righty_vert[righty_idx[i]]));
        }
        SortSync(righty_betweens, righty_idx, righty_vals);
        for ( size_t i = 0; i < righty_idx.size(); ++i ) {
            cout << "righty " << i << " betweens=" << righty_betweens[i] <<
                    ", vert=" << graph.Vert(righty_vert[righty_idx[i]]) << ", vals=" <<
                    righty_vals[i] << endl;
        }
        // break ties
        int j = 1;
        while ( j < righty_betweens.isize() && righty_betweens[j] == righty_betweens[0] ) j++;
        if ( j < righty_betweens.isize() ) {
            righty_vals.erase( righty_vals.begin()+j, righty_vals.end() );
            righty_idx.erase( righty_idx.begin()+j, righty_idx.end() );
        }
        SortSync(righty_vals, righty_idx );
        for ( size_t i = 0; i < righty_idx.size(); ++i ) {
            cout << "righty " << i << " betweens=" << righty_betweens[i] <<
                    ", vert=" << graph.Vert(righty_vert[righty_idx[i]]) << ", vals=" <<
                    righty_vals[i] << endl;
        }

        right = graph.Vert(righty_vert[righty_idx[0]]);
        found_right = true;
    }


    return make_pair(found_left, found_right );
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
    CommandArgument_IntSet_OrDefault_Doc(TRACE, "",
            "list of lines to trace in some way");
    CommandArgument_Double_OrDefault_Doc(VTHRESH, 1e-5,
            "ratio threshold above which edges will not be output in final graph.");
    CommandArgument_Int(DEPTH);
    CommandArgument_Bool_OrDefault(VERBOSE,True);
    CommandArgument_String_OrDefault_Doc(DUMPHEAD,"",
            "dump some debugging files with this head");
    CommandArgument_String_OrDefault_Doc(OUTHEAD,"",
            "output graph head" );
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
    testLines.push_back(cseed);
    auto rateThresh = rates[0].mHitsPerKb * 0.10;
    for ( auto const& hit : rates ) {
        auto lineId = hit.mLineId;
        auto lineLen = hic.mLineInfoVec[lineId].mLineLen;
        if ( lineLen < 6000 ) continue;
//        if ( testLines.isize() == DEPTH || hit.mHitsPerKb < rateThresh ) break;
        if ( testLines.isize() == DEPTH+1 ) break;      // +1 because the seed is there
        cout << hit << ", ";
        testLines.push_back( lineId );
    }
    cout << endl;


#if 0
    // checking whether cohort is canonical
    for ( auto const& test : testLines )
        cout << test << " => " << canonicalLine(test,hic.mLineInfoVec) <<
        "|~" << hic.mLineInfoVec[test].mLineIdRC << endl;
#endif // ...it is

    // pull hiCPairs incident on these lines
    cout << Date() << ": reading HiC pairs for selected lines" << endl;
    HICVec hiCPairs;
    getHiCPairs( HICPAIRS, hic.mEdgeToLineMap, &hiCPairs,
             false /* don't skip intra-line pairs */, testLines );

    cout << Date() << ": done" << endl;

    // count hicPairs by line if trace
    if ( TRACE.size() ) {
        size_t count = 0;
#pragma omp parallel for reduction(+:count)
        for ( size_t i = 0; i < hiCPairs.size(); ++i ) {
            auto const& hiCPair = hiCPairs[i];
            int line1 = hic.mEdgeToLineMap[hiCPair.el1().getEdgeId()];
            int line2 = hic.mEdgeToLineMap[hiCPair.el2().getEdgeId()];
            if ( Member(TRACE, line1) && Member(TRACE, line2) )
                count++;
        }
        cout << Date() << ": TRACE on hiCPairs count=" << count << endl;
    }

    cout << Date() << ": calculating edge pairs -> line pairs" << endl;
    LineOffsetCalc offCalc( hic.mHBX, hic.mLines, hic.mInv,
            hic.mEdgeToLineMap, hic.mRuler, LineOffsetCalc::THREAD_SAFE );

    cout << Date() << ": preloading cache of line lengths and edge<->line map" << endl;
    for ( auto const& testLine : testLines )
        offCalc.cachePreload(testLine);

    LineLocPairVec linePairs(hiCPairs.size());
#pragma omp for
    for ( size_t i = 0; i < hiCPairs.size(); ++i ) {
        auto const & hiCPair = hiCPairs[i];
        int line1 = hic.mEdgeToLineMap[hiCPair.el1().getEdgeId()];
        int line2 = hic.mEdgeToLineMap[hiCPair.el2().getEdgeId()];
        if ( !Member(testLines, line1 ) ||
                !Member(testLines, line2 ) ) continue;
        LineLocPair tmp( offCalc(hiCPair.el1()), offCalc(hiCPair.el2()) );
        linePairs[i]=tmp;
    }

    if ( DUMPHEAD != "" )
        BinaryWriter::writeFile(DUMPHEAD+".llpairs", linePairs);

    cout << "linePairs.size()=" << linePairs.size() << endl;

    if ( TRACE.size() ) {
        size_t count = 0;
#pragma omp parallel for reduction(+:count)
        for ( size_t i = 0; i < linePairs.size(); ++i ) {
            auto const& linePair = linePairs[i];
            if ( Member(TRACE, linePair.ll1().getLineId()) &&
                    Member(TRACE, linePair.ll2().getLineId() ) )
                count++;
        }
        cout << Date() << ": TRACE on linePairs count=" << count << endl;
    }

    // length of the lines in question
    if ( VERBOSE ) { cout << "Line Lengths: " << endl; }
    std::unordered_map<int, int> testLineLengths;
    for ( auto const lineId : testLines ) {
        testLineLengths[lineId] = GetLineLength(hic.mLines[lineId], hic.mRuler);
        if ( VERBOSE) {cout << "lineId " << lineId << ": " << testLineLengths[lineId] << endl;}
    }

    // calculate all pairwise order and orientations
    digraphVE<TaggedLineId,EdgeEvidence> graph;
    vec<String> vertex_labels;
    vec<String> vertex_colors;
    for ( size_t i = 0; i < testLines.size(); ++i ) {
        graph.AddVertex(TaggedLineId(testLines[i], false));     // A
        vertex_labels.push_back(ToString(testLines[i])+" ");
        vertex_colors.push_back((i==0) ? "green":"black");
        graph.AddVertex(TaggedLineId(testLines[i], true));      // A'
        vertex_labels.push_back(ToString(testLines[i])+"'");
        vertex_colors.push_back((i==0) ? "green":"black");
    }

    vec<vec<String>> edge_labels(graph.N());

    int nLines = testLines.isize();
#pragma omp parallel for if(!VERBOSE)
    for ( int z = 0; z < nLines*nLines; z++ ) {
        int i = z / nLines;
        int j = z % nLines;
        if ( j <= i ) continue;

        auto const& A = testLines[i];
        auto const& B = testLines[j];

        auto scores = testLinePair(TaggedLineId(A,false), TaggedLineId(B,false), testLineLengths, linePairs, hic.mSepDist);
        if ( VERBOSE ) {
            cout << "A=[" << A << "], B=[" << B << "] => ";
        }
        auto win = winningPairIndex(scores);
        if ( scores[win.first].mCount == 0 ) {
            if ( VERBOSE ) {
                cout << "no evidence" << endl;
            }
            continue;
        }
        auto ratio = winningPairRatio(scores, Scores::LOG_LIK_A);
        String config = pairIndexToConfigString(win.first);
        String config_alt = pairIndexToConfigString(win.second);
        if ( VERBOSE ) {
            cout << config << " ratio=" << ratio << ", alt=" << config_alt << ", count=" << scores[0].mCount << endl;
        }
        auto pair = pairIndexToConfig(win.first, TaggedLineId(A,false), TaggedLineId(B,false));

        int v, w;
        if ( pair.first.getId() == A && pair.second.getId() == B ) {
            v = 2*i;
            w = 2*j;
        } else if ( pair.first.getId() == B && pair.second.getId() == A ) {
            v = 2*j;
            w = 2*i;
        } else FatalErr("BUG");

        if ( pair.first.isTag() ) v++;
        if ( pair.second.isTag() ) w++;

#pragma omp critical
        graph.AddEdge(v,w, EdgeEvidence{ratio,scores[win.first].mCount,win.second});
    }

    // okay, infer left and right adjacencies based on graph info
    TaggedLineId left, right;
    auto indicator = inferLeftRight2( cseed, graph, left, right );
    cout << "final config = " << "[";
    if ( indicator.first )
        cout << left.getId()
            << ( left.isTag() ? "'" : "" );
    else
        cout << "X";
    cout << "] " << cseed << " [";
    if ( indicator.second )
        cout << right.getId()
        << (right.isTag() ? "'" : "" );
    else
        cout << "X";
    cout << "]" << endl;

    if ( OUTHEAD != "" ) {
        // for visualization purposes, delete the non-passing edges
        for ( int v = 0; v < graph.N(); ++v ) {
            vec<int> to_delete;
            for ( int j = 0; j < graph.FromSize(v); ++j )
                if ( graph.EdgeObjectByIndexFrom(v,j).ratio > VTHRESH ) {
                    to_delete.push_back(j);
                }
            graph.DeleteEdgesFrom(v, to_delete);
            EraseTheseIndices(edge_labels[v],to_delete);
        }

        // don't print singletons
        for ( int v = 0; v < graph.N(); ++v ) {
            if ( graph.FromSize(v) == 0 && graph.ToSize(v) == 0 )
                vertex_labels[v] = "";
        }

        dumpGraph(graph,"graph.dot", cseed );
    }

    return 0;
}


#include "graph/DigraphTemplate.h"
template void digraphVE<TaggedLineId,EdgeEvidence>::AddVertex(TaggedLineId const&);
template TaggedLineId const& digraphVE<TaggedLineId,EdgeEvidence>::Vert(int) const;

