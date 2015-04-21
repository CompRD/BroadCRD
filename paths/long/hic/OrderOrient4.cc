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
#include "paths/long/hic/Scores.h"
#include "graph/Digraph.h"
#include "paths/long/hic/LineGraph.h"
#include "paths/long/hic/CanonicalTestOrder.h"
#include "paths/long/hic/LineMapReader.h"
#include "feudal/BinaryStream.h"
#include "paths/long/hic/ScaffoldPile.h"
#include "paths/long/hic/NeighborScaff.h"
#include <array>


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

#if 0
double winningRatio(vec<Scores> scores, Scores::Metric metric)
{
    auto idx = winningPairIndex(scores).first;

#if 0
    if ( metric == Scores::LOG_LIK_A || metric == Scores::LOG_LIK ) {
       auto count =  scores[0].mCount;
       if ( idx == 0 || idx == 2 ) return exp( ab - ba ) * count / 100000;
       else return exp( ba - ab ) * count / 100000;
    } else {
        if ( idx == 0 || idx == 2 ) return (ba/ab); return (ab/ba);
    }
#endif
}
#endif

Bool acyclic( LineGraph const& graph0,
        int cseed,
        const double upper_thresh = std::numeric_limits<double>::max() )
{
    LineGraph graph(graph0);

    for ( int v = 0; v < graph.N(); ++v ) {
        vec<int> to_delete;
        for ( int j = 0; j < graph.FromSize(v); ++j )
            if ( graph.EdgeObjectByIndexFrom(v,j).ratio > upper_thresh ) {
                to_delete.push_back(j);
            }
        graph.DeleteEdgesFrom(v, to_delete);
    }

    // remove components not containing the seed
    vec<vec<int>> components;
    graph.Components(components);
    if ( components.size() > 1 ) {
        vec<int> to_delete;
        for ( size_t i = 0; i < components.size(); ++i ) {
            bool found = false;
            for ( auto v : components[i] ) {
                if ( graph.Vert(v) == TaggedLineId( cseed, false ) ) {
                    found = true;
                    break;
                }
            }
            if ( !found ) to_delete.append(components[i]);
        }
        for ( auto v : to_delete )
            graph.DeleteEdgesAtVertex(v);
    }

    return graph.Acyclic();
}

size_t pathsBetweenWithout( LineGraph const& graph, int v, int w, int without_edge )
{
    ForceAssertNe(v,w);

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
            for ( int j = 0; j < graph.FromSize(v); ++j ) {
                auto const& edge = graph.EdgeObjectByIndexFrom(v,j);
                int edgei = graph.EdgeObjectIndexByIndexFrom(v,j);
                int from = graph.From(v)[j];
                if ( edgei != without_edge )
                    stack.push_back( from );
            }
        }
    }

    return hits;
}

size_t pathsBetween( LineGraph const& graph, int v, int w,
        const double upper_thresh = std::numeric_limits<double>::max() )
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
            for ( int j = 0; j < graph.FromSize(v); ++j ) {
                auto const& edge = graph.EdgeObjectByIndexFrom(v,j);
                int from = graph.From(v)[j];
                if ( edge.ratio <= upper_thresh )
                    stack.push_back( from );
            }
        }
    }

    return hits;
}

double maximalAcyclicRatio(LineGraph const& graph, int cseed, double start,
        double end, double step = 1)
{
    for ( double thresh = start; thresh >= end; thresh -= step ) {
        if (acyclic(graph, cseed, thresh)) return thresh;
    }

    return start;
}


std::pair<double,double> inferLeftRight2( HiCExperiment const& hic,
        int cseed, LineGraph const& graph0,
        TaggedLineId& left, TaggedLineId& right, int iteration, Bool VERBOSE = False )
{
    double const start_thresh=1e-25;
    double const end_thresh=1e-10;
    double const delta_thresh=10;


    ForceAssertGt(graph0.N(),0);
    auto seed = graph0.Vert(0);

    double left_ratio = 0., right_ratio = 0.;

    LineGraph graph(graph0);

    // freaky test - remove all sources and sinks
    vec<int> to_delete;
#if 1
#warning "freaky test - remove all sources and sinks on a single stalk"
    for ( int v = 0; v < graph.N(); ++v )
        if ( graph.FromSize(v)+graph.ToSize(v) <= 1) to_delete.push_back(v);
    for ( auto v : to_delete )
        graph.DeleteEdgesAtVertex(v);
#endif


    // first find the left threshold
#if 0
    double thresh = start_thresh;
    for ( Bool done = False ; !done && thresh <= end_thresh; thresh *= delta_thresh ) {
        for ( int j = 0; j < graph.ToSize(0); ++j ) {
            size_t npaths = pathsBetween( graph, graph.To(0)[j], 0, thresh );
//            PRINT3(npaths, thresh, graph.Vert(graph.To(0)[j]).getId() );
            if ( npaths > 0 ) {
                done = True;
                break;
            }
        }
    }
#endif
    double tstart = -10.;
    double tend = -55.;

    double thresh = maximalAcyclicRatio(graph,cseed,tstart,tend);
//    if ( thresh == 0. ) return std::make_pair(false,false);
//    PRINT(thresh);

    if ( VERBOSE ) cout << "*** THRESHOLD CHOSEN = " << thresh << endl;

    // list left nodes and confidences
    vec<int> lefty_vidx;        // index only into lefty_vert
    vec<double> lefty_vals;
    vec<int> lefty_vert;
//    double thresh;
//    for ( thresh = start_thresh;
//            thresh <= end_thresh; thresh *= delta_thresh ) {
        for ( int j = 0; j < graph.ToSize(0); ++j ) {
            auto ratio = graph.EdgeObjectByIndexTo(0,j).ratio;
//            double len = hic.mLineInfoVec[graph.Vert(graph.To(0)[j]).getId()].mLineLen;
#warning "experimental"
            double len = graph.EdgeObjectByIndexTo(0,j).sumLen;
            double val = graph.EdgeObjectByIndexTo(0,j).count / len;
            if ( ratio < thresh ) {
#warning "also experimental"
                lefty_vals.push_back( ratio );
//                lefty_vals.push_back(-val);
                lefty_vidx.push_back(j);
                if (VERBOSE) cout << "accepting lefty " << graph.Vert(graph.To(0)[j]) << ", ratio=" << ratio << endl;
            }
            else {
                if (VERBOSE) cout << "failing lefty " << graph.Vert(graph.To(0)[j]) << ", ratio=" << ratio << endl;
            }
            lefty_vert.push_back(graph.To(0)[j]);
        }
//        if ( lefty_vals.size() > 0 ) break;
//        else lefty_vert.clear();
//    }

//      PRINT(lefty_idx.size());

    // find winner on the left
#if 0   // original scoring: betweens -> tie breaking with value
    if ( lefty_idx.size() > 0 ) {
        vec<int> lefty_betweens;
        for ( size_t i = 0; i < lefty_idx.size(); ++i ) {
            lefty_betweens.push_back( pathsBetween(graph, lefty_vert[lefty_idx[i]], 0, thresh) );
        }
        SortSync(lefty_betweens, lefty_idx, lefty_vals);
        for ( size_t i = 0; i < lefty_idx.size(); ++i ) {
            if ( VERBOSE ) cout << "lefty " << i << " betweens=" << lefty_betweens[i] <<
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
            if (VERBOSE) cout << "lefty " << i << " betweens=" << lefty_betweens[i] <<
                    ", vert=" << graph.Vert(lefty_vert[lefty_idx[i]]) << ", vals=" <<
                    lefty_vals[i] << endl;
        }

        left = graph.Vert(lefty_vert[lefty_idx[0]]);
        found_left = true;
    }

#else
    if ( lefty_vidx.size() > 0 ) {
        vec<int> lefty_betweens;
        for ( size_t i = 0; i < lefty_vidx.size(); ++i ) {
            lefty_betweens.push_back( pathsBetween(graph, lefty_vert[lefty_vidx[i]], 0, thresh) );
        }
        vec<double> lefty_vals2(lefty_vals);
        for ( size_t i = 0; i < lefty_vals2.size(); ++i ) {
            lefty_vals2[i] = lefty_vals2[i] - ( lefty_betweens[i] - 1 ) * lefty_vals2[i] * 0.40;
        }
        SortSync(lefty_vals2, lefty_vidx, lefty_vals, lefty_betweens);
        for ( size_t i = 0; i < lefty_vidx.size(); ++i ) {
            if (VERBOSE) cout << "lefty " << i << " betweens=" << lefty_betweens[i] <<
                    ", vert=" << graph.Vert(lefty_vert[lefty_vidx[i]]) << ", vals=" <<
                    lefty_vals[i] << ", vals2=" << lefty_vals2[i] << endl;
        }
        // less than because it's negative in the next line
        if ( lefty_vidx.size() == 1  || lefty_vals2[0] < lefty_vals2[1] * 1.05 ) {
            left = graph.Vert(lefty_vert[lefty_vidx[0]]);
            left_ratio = lefty_vals[0];
        }
    }
#endif


    // now find the right threshold
#if 0
    thresh = start_thresh;
    for ( Bool done = False ; !done && thresh <= end_thresh; thresh *= delta_thresh ) {
        for ( int j = 0; j < graph.ToSize(0); ++j ) {
            size_t npaths = pathsBetween( graph, 0, graph.From(0)[j], thresh );
//            PRINT3(npaths, thresh, graph.Vert(graph.To(0)[j]).getId() );
            if ( npaths > 0 ) {
                done = True;
                break;
            }
        }
    }
#endif
    // list right nodes and confidences
    vec<int> righty_vidx;       // index only into righty_vert
    vec<double> righty_vals;
    vec<int> righty_vert;
    /* double thresh -- common with above */
//    for ( thresh = start_thresh;
//            thresh <= end_thresh; thresh *= delta_thresh ) {
        for ( int j = 0; j < graph.FromSize(0); ++j ) {
            auto ratio = graph.EdgeObjectByIndexFrom(0,j).ratio;
//            double len = hic.mLineInfoVec[graph.Vert(graph.From(0)[j]).getId()].mLineLen;
#warning "experimental"
            double len = graph.EdgeObjectByIndexFrom(0,j).sumLen;
            double val = graph.EdgeObjectByIndexFrom(0,j).count / len;
            if ( ratio < thresh ) {
#warning "also experimental"
                righty_vals.push_back( ratio );
//                righty_vals.push_back( -val );
                righty_vidx.push_back(j);
                if (VERBOSE) cout << "accepting righty " << graph.Vert(graph.From(0)[j]) << ", ratio=" << ratio << endl;
            } else {
                if (VERBOSE) cout << "failing righty " << graph.Vert(graph.From(0)[j]) << ", ratio=" << ratio << endl;
            }
            righty_vert.push_back(graph.From(0)[j]);
        }
//        if ( righty_vals.size() > 0 ) break;
//        else righty_vert.clear();
//    }

    // find winner on the right
#if 0
    if ( righty_idx.size() > 0 ) {
        vec<int> righty_betweens;
        for ( size_t i = 0; i < righty_idx.size(); ++i ) {
            righty_betweens.push_back( pathsBetween(graph,0,righty_vert[righty_idx[i]], thresh));
        }
        SortSync(righty_betweens, righty_idx, righty_vals);
        for ( size_t i = 0; i < righty_idx.size(); ++i ) {
            if ( VERBOSE ) cout << "righty " << i << " betweens=" << righty_betweens[i] <<
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
            if (VERBOSE) cout << "righty " << i << " betweens=" << righty_betweens[i] <<
                    ", vert=" << graph.Vert(righty_vert[righty_idx[i]]) << ", vals=" <<
                    righty_vals[i] << endl;
        }

        right = graph.Vert(righty_vert[righty_idx[0]]);
        found_right = true;
    }
#else
    if ( righty_vidx.size() > 0 ) {
        vec<int> righty_betweens;
        for ( size_t i = 0; i < righty_vidx.size(); ++i ) {
            righty_betweens.push_back( pathsBetween(graph, 0, righty_vert[righty_vidx[i]], thresh) );
        }
        vec<double> righty_vals2(righty_vals);
        for ( size_t i = 0; i < righty_vals2.size(); ++i ) {
            righty_vals2[i] = righty_vals2[i] - ( righty_betweens[i] - 1 ) * righty_vals2[i] * 0.40;
        }
        SortSync(righty_vals2, righty_vidx, righty_vals, righty_betweens);
        for ( size_t i = 0; i < righty_vidx.size(); ++i ) {
            if (VERBOSE) cout << "righty " << i << " betweens=" << righty_betweens[i] <<
                    ", vert=" << graph.Vert(righty_vert[righty_vidx[i]]) << ", vals=" <<
                    righty_vals[i] << ", vals2=" << righty_vals2[i] << endl;
        }
        // less than because it's negative in the next line
        if ( righty_vidx.size() == 1 || righty_vals2[0] < righty_vals2[1] * 1.05 ) {
            right = graph.Vert(righty_vert[righty_vidx[0]]);
            right_ratio = righty_vals[0];
        }
    }
#endif


    return make_pair(left_ratio, right_ratio);
}


int neighborLines( int seed, HiCExperiment const& hic,
        int DEPTH,
        vec<int>& testLines,
        Bool verb = False)
{
    auto cseed = canonicalLine(seed, hic.mLineInfoVec);
    if ( verb) cout << "HitRateVec[" << cseed << "]=";
    auto rates = hic.mHitRates[cseed];
    if ( rates.size() == 0 ) {
        if (verb) cout << "nothing to do!" << endl;
        return cseed;
    }
    testLines.push_back(cseed);
    auto rateThresh = rates[0].mHitsPerKb * 0.10;
    for ( auto const& hit : rates ) {
        auto lineId = hit.mLineId;
        auto lineLen = hic.mLineInfoVec[lineId].mLineLen;
        if ( lineLen < 6000 ) continue;
//        if ( testLines.isize() == DEPTH || hit.mHitsPerKb < rateThresh ) break;
        if ( testLines.isize() == DEPTH+1 ) break;      // +1 because the seed is there
        if (verb) cout << hit << ", ";
        testLines.push_back( lineId );
    }
    if (verb) cout << endl;
    return cseed;
}


template <class V>
void flipScaff(V& scaff, vec<int> lineInv)
{
    for ( auto& line : scaff ) {
        ForceAssertNe(lineInv[line], -1);
        line = lineInv[line];
    }
    std::reverse(scaff.begin(), scaff.end());
}

void flipTaggedArr(std::array<TaggedLineId,3>& scaff)
{
    for ( auto& line : scaff )
        if (line.getId() != TaggedLineId::maxVal ) line.flipTag();
    std::swap(scaff[0], scaff[2]);
}

void flipTaggedLineVec(vec<TaggedLineId>& scaff)
{
    for ( auto& line : scaff ) line.flipTag();
    std::reverse(scaff.begin(), scaff.end());
}

size_t scaffLen( vec<int> scaff, LineInfoVec const& linfo )
{
    size_t len = std::accumulate( scaff.begin(), scaff.end(), 0u,
            [&linfo](size_t len, int lineid) {
                return len + linfo[lineid].mLineLen; });
    return len;
}

std::array<TaggedLineId,3> hiCGrowNeighborhoodScaff(int seed,
        ScaffoldPile const& scaffPile,
        LineLocPairVec const& linePairs0,
        vec<int> const& lineInv,
        unordered_map<int,int> const& lineLengths,
        HiCExperiment const& hic,
        int DEPTH,
        int MINLINESIZE,
        int iter,
        Bool VERBOSE = False,
        String OUTHEAD = "",
        double VTHRESH = 1e-5 )
{
    const size_t HEAD = 2;      // how many lines at ends of each scaffold to use for seeding

    std::array<TaggedLineId,3> triplet;
    triplet[0] = triplet[2] = TaggedLineId(TaggedLineId::maxVal,false);
//    triplet[1] = TaggedLineId(seed, false);

    // turn seed into scaffold start
    int start_seed = scaffPile.start(seed);
    triplet[1] = TaggedLineId(start_seed, false);
    if ( VERBOSE ) cout << Date() << ": calculating cohort based on seed "
            << start_seed << ": " << endl;

    // figure out which scaffolds are adjacent
    vec<int> testStarts;
    bool verb = false;
//    if ( start_seed == 6841 ) verb = true;
    neighborScaff( start_seed, scaffPile, hic, DEPTH, MINLINESIZE, testStarts, HEAD, verb );
    if ( VERBOSE ) cout << Date() << printSeq(testStarts) << endl;

    // turn testStarts into testLines w/ inverses so we can pick out line pairs
    vec<int> testLines;
    for ( auto const start : testStarts ) {
        vec<int> scaff = scaffPile(start);
        for ( int line : scaff ) {
            testLines.push_back(line);
            int invLine = lineInv[line];
            testLines.push_back(invLine);
        }
    }

//    cout << "neighborLines for seed(start) " << seed << "(" << start_seed
//            << ")=" << printSeq(testLines) << endl;

    // pull out pairs from lines in question
    LineLocPairVec linePairs;
#pragma omp parallel for
    for ( size_t i = 0; i < linePairs0.size(); ++i ) {
        auto const& pair  = linePairs0[i];
        if ( Member(testLines, pair.first.getLineId()) &&
                Member(testLines, pair.second.getLineId()) )
#pragma omp critical
            linePairs.push_back(pair);
    }

    // calculate all pairwise order and orientations
    digraphVE<TaggedLineId,EdgeEvidence> graph;
    for ( size_t i = 0; i < testStarts.size(); ++i ) {
        graph.AddVertex(TaggedLineId(testStarts[i], false));     // A
        graph.AddVertex(TaggedLineId(testStarts[i], true));      // A'
    }

    int nScaffs = testStarts.isize();
    size_t nthreads = omp_get_num_threads();
//    if ( VERBOSE ) omp_set_num_threads(1);
//#pragma omp parallel for
    for ( int z = 0; z < nScaffs*nScaffs; z++ ) {
        std::ostringstream out;
        int i = z / nScaffs;
        int j = z % nScaffs;
        if ( j <= i ) continue;

        int A = testStarts[i];
        int B = testStarts[j];

        int Alen = scaffLen( scaffPile(A), hic.mLineInfoVec );
        int Blen = scaffLen( scaffPile(B), hic.mLineInfoVec );

        auto scores = testScaffPair(scaffPile, A, B, lineLengths, linePairs, hic.mSepDist, hic.mLineInfoVec);
        if ( VERBOSE ) {
            out << "A=[" << printSeq(scaffPile(A)) << "], B=[" << printSeq(scaffPile(B)) << "] => ";
        }
//        auto win_old = winningPairIndex(scores);
        double ratio;
        auto win = winningPairOrderOrient(scores,ratio);
//        auto win = winningPairScore(scores,ratio);
//        PRINT(ratio);
        if ( scores[win.first].mCount == 0 || ratio == 0. ) {
            if ( VERBOSE ) {
                out << "no evidence" << endl;
            }
            continue;
        }
        Bool verb=VERBOSE;
        if ( verb ) {
            out << "Testing " << printSeq(scaffPile(A)) << " versus " << printSeq(scaffPile(B)) << endl;
        }
//        auto ratio_old = winningOrderRatio(scores, Scores::LOG_LIK_A, verb);
//        PRINT2(win.first,win.second);
//        PRINT2(win_old.first,win_old.second);
//        PRINT2(ratio, ratio_old);
        String config = CanonicalTestOrder::pairIndexToConfigString(win.first);
        String config_alt = CanonicalTestOrder::pairIndexToConfigString(win.second);
        if ( VERBOSE ) {
            out << config << " ratio=" << ratio << ", alt=" << config_alt << ", count=" << scores[0].mCount
                    << ", count/len=" << scores[0].mCount / static_cast<double>(Alen+Blen) << endl;
            for ( size_t i = 0; i < scores.size(); ++i )
                out << "\t" << CanonicalTestOrder::pairIndexToConfigString(i) << "=>" <<
                        scores[i](Scores::LOG_LIK_A) << " | ";
            out << endl;
        }
        auto pair = CanonicalTestOrder::pairIndexToConfig(win.first, TaggedLineId(A,false), TaggedLineId(B,false));

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

        size_t len = Alen+Blen;

#pragma omp critical
        {
        graph.AddEdge(v,w, EdgeEvidence{ratio,scores[win.first].mCount, len, win.second});
        cout << out.str();
        }
    }

    // okay, infer left and right adjacencies based on graph info
    TaggedLineId left, right;
    auto ratios = inferLeftRight2( hic, seed, graph, left, right, iter, VERBOSE );

    double thresh = 0.;

    if ( ratios.first != 0. && ratios.first < thresh ) triplet[0] = left;
    if ( ratios.second != 0. && ratios.second < thresh ) triplet[2] = right;

//    if ( cseed != seed ) flipScaff(triplet, lineInv);


    if ( VERBOSE ) {
        cout << "final config = " << "[";
        if ( ratios.first != 0. )
            cout << left.getId()
            << ( left.isTag() ? "'" : "" );
        else
            cout << "X";

        cout << "] " << seed << " [";
        if ( ratios.second != 0. )
            cout << right.getId()
            << (right.isTag() ? "'" : "" );
        else
            cout << "X";
        cout << "]" << endl;
    }


    if ( OUTHEAD != "" ) {
        // for visualization purposes, delete the non-passing edges
        for ( int v = 0; v < graph.N(); ++v ) {
            vec<int> to_delete;
            for ( int j = 0; j < graph.FromSize(v); ++j )
                if ( graph.EdgeObjectByIndexFrom(v,j).ratio > VTHRESH ) {
                    to_delete.push_back(j);
                }
            graph.DeleteEdgesFrom(v, to_delete);
        }

        dumpLineGraph(graph,OUTHEAD+"."+ToString(seed)+"."+ToString(iter)+".dot", seed );
    }

    return triplet;
}

std::array<int,3> hiCGrowNeighborhood(int seed,
        LineLocPairVec const& linePairs0,
        vec<int> const& lineInv,
        HiCExperiment const& hic,
        int DEPTH,
/*        TaggedLineId& left, TaggedLineId& right, */
        Bool VERBOSE = False,
        String OUTHEAD = "",
        double VTHRESH = 1e-5, ScaffoldPile* scaffp = nullptr )
{
    std::array<int,3> triplet;
    triplet[0] = triplet[2] = -1; triplet[1] = seed;
    if ( VERBOSE ) cout << Date() << ": calculating cohort based on seed" << endl;
    vec<int> testLines;
    int cseed = neighborLines( seed, hic, DEPTH, testLines );

    if ( scaffp )
    {
        vec<int> starts;
        neighborScaff( cseed, *scaffp, hic, DEPTH, 6000, starts);
#if 0
        vec<int> testLines2(testLines);
        Sort(starts);
        Sort(testLines2);
        if ( starts != testLines2 ) {
            cout << "mismatch for seed " << seed << endl;
            cout << "starts: " << printSeq(starts) << endl;
            cout << "testLines2: " << printSeq(testLines2) << endl;
        } else
            cout << "ALL IS WELL" << endl;
#endif
        testLines = starts;
    }


    if ( testLines.size() == 0 ) return triplet;

    if ( VERBOSE ) {
        if ( cseed != seed ) cout <<
                "INFO: seed canonicalized to " << cseed << endl;
        cout << "neighbors[" << cseed << "]=";
        cout << printSeq(testLines) << endl;
    }

    // length of the lines in question
    std::unordered_map<int, int> testLineLengths;
    for ( auto const lineId : testLines ) {
        testLineLengths[lineId] = hic.mLineInfoVec[lineId].mLineLen;

#ifdef EXTRATESTING
        auto tmp = GetLineLength(hic.mLines[lineId], hic.mRuler);
        ForceAssertEq(tmp,testLineLengths[lineId]);
#endif
    }

    // pull out pairs from lines in question
//    cout << Date() << ": pulling out line pairs" << endl;
    LineLocPairVec linePairs;
#pragma omp parallel for
    for ( size_t i = 0; i < linePairs0.size(); ++i ) {
        auto const& pair  = linePairs0[i];
        if ( Member(testLines, pair.first.getLineId()) &&
                Member(testLines, pair.second.getLineId()) )
#pragma omp critical
            linePairs.push_back(pair);
    }
//    cout << Date() << ": done" << endl;
#if 0
    for ( auto const testLine : testLines )
        for ( auto const pairId : lpIndex[testLine] ) {
            auto pair = linePairs0[pairId];
            if ( Member(testLines, pair.first.getLineId()) &&
                    Member(testLines, pair.second.getLineId()) )
                linePairs.push_back( linePairs0[pairId] );
        }
#endif
//    PRINT(linePairs.size());

    // calculate all pairwise order and orientations
    digraphVE<TaggedLineId,EdgeEvidence> graph;
    for ( size_t i = 0; i < testLines.size(); ++i ) {
        graph.AddVertex(TaggedLineId(testLines[i], false));     // A
        graph.AddVertex(TaggedLineId(testLines[i], true));      // A'
    }

    int nLines = testLines.isize();
    size_t nthreads = omp_get_num_threads();
//    if ( VERBOSE ) omp_set_num_threads(1);
#pragma omp parallel for
    for ( int z = 0; z < nLines*nLines; z++ ) {
        std::ostringstream out;
        int i = z / nLines;
        int j = z % nLines;
        if ( j <= i ) continue;

        auto const& A = testLines[i];
        auto const& B = testLines[j];

        int Alen = hic.mLineInfoVec[A].mLineLen;
        int Blen = hic.mLineInfoVec[B].mLineLen;

        auto scores = testLinePair(TaggedLineId(A,false), TaggedLineId(B,false), testLineLengths, linePairs, hic.mSepDist);
        if ( VERBOSE ) {
            out << "A=[" << A << "], B=[" << B << "] => ";
        }
        auto win = winningPairIndex(scores);
        if ( scores[win.first].mCount == 0 ) {
            if ( VERBOSE ) {
                out << "no evidence" << endl;
            }
            continue;
        }
//        auto ratio = winningPairRatio(scores, Alen+Blen, Scores::LOG_LIK_A);
        auto ratio = winningOrderRatio(scores, Scores::LOG_LIK_A);
        String config = CanonicalTestOrder::pairIndexToConfigString(win.first);
        String config_alt = CanonicalTestOrder::pairIndexToConfigString(win.second);
        if ( VERBOSE ) {
            out << config << " ratio=" << ratio << ", alt=" << config_alt << ", count=" << scores[0].mCount
                    << ", count/len=" << scores[0].mCount / static_cast<double>(Alen+Blen) << endl;
            for ( size_t i = 0; i < scores.size(); ++i )
                out << "\t" << CanonicalTestOrder::pairIndexToConfigString(i) << "=>" <<
                        scores[i](Scores::LOG_LIK_A) << " | ";
            out << endl;
        }
        auto pair = CanonicalTestOrder::pairIndexToConfig(win.first, TaggedLineId(A,false), TaggedLineId(B,false));

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

        size_t len = Alen+Blen;

#pragma omp critical
        {
        graph.AddEdge(v,w, EdgeEvidence{ratio,scores[win.first].mCount, len, win.second});
        cout << out.str();
        }
    }

    // okay, infer left and right adjacencies based on graph info
    TaggedLineId left, right;
    auto ratios = inferLeftRight2( hic, cseed, graph, left, right, -1, VERBOSE );

    // handle any RC edges
    if ( ratios.first != 0. ) triplet[0] =
            left.isTag() ? lineInv[left.getId()] : left.getId();
    triplet[1] = cseed;
    if ( ratios.second != 0. ) triplet[2] =
            right.isTag() ? lineInv[right.getId()] : right.getId();

    if ( cseed != seed ) flipScaff(triplet, lineInv);


    if ( VERBOSE ) {
        cout << "final config = " << "[";
        if ( ratios.first != 0. )
            cout << left.getId()
            << ( left.isTag() ? "'" : "" );
        else
            cout << "X";
        cout << "] " << cseed << " [";
        if ( ratios.second != 0. )
            cout << right.getId()
            << (right.isTag() ? "'" : "" );
        else
            cout << "X";
        cout << "]" << endl;
    }


    if ( OUTHEAD != "" ) {
        // for visualization purposes, delete the non-passing edges
        for ( int v = 0; v < graph.N(); ++v ) {
            vec<int> to_delete;
            for ( int j = 0; j < graph.FromSize(v); ++j )
                if ( graph.EdgeObjectByIndexFrom(v,j).ratio > VTHRESH ) {
                    to_delete.push_back(j);
                }
            graph.DeleteEdgesFrom(v, to_delete);
        }

        dumpLineGraph(graph,OUTHEAD+"."+ToString(cseed)+".dot", cseed );
    }

    return triplet;
}

template < typename T >
void printCols( size_t ncols, T const& s )
{
    std::ostringstream out;
    out << s;
    while ( out.str().size() < ncols ) {
        out << " ";
    }
    cout << out.str();
}


void addToLineGraph( TaggedLineId tv, TaggedLineId tw, LineGraph& graph,
        std::map<TaggedLineId, int>& mapLineToVertex )
{
    int v;
    try { v = mapLineToVertex.at( tv ); }
    catch ( std::out_of_range const& e ) {
        v = graph.N();
        graph.AddVertex(tv);
        mapLineToVertex[tv] = v;
    }
    int w;
    try { w = mapLineToVertex.at( tw ); }
    catch ( std::out_of_range const& e ) {
        w = graph.N();
        graph.AddVertex(tw);
        mapLineToVertex[tw] = w;
    }
    graph.AddEdge(v,w,EdgeEvidence());
}


void testEndsForCompatibility( TaggedLineId end1, TaggedLineId end2,
        LineGraph const& graph, HiCExperiment const& hic,
        std::map<TaggedInt, int> const& mapLineToVertex,
        Bool VERBOSE )
{
    const int backup = 5;

    int v1 = mapLineToVertex.at(end1);
    int v2 = mapLineToVertex.at(end2);

    bool v1rc = false, v2rc = false;

}



void cleanAmbiguities( LineGraph& graph, HiCExperiment const& hic, Bool VERBOSE )
{
    // dumb dumb dumb transitive reduction
    for ( int v = 0; v < graph.N(); ++v ) {
        vec<int> to_delete;
        for ( int j = 0; j < graph.FromSize(v); ++j ) {
            int w = graph.From(v)[j];
            int edgei = graph.EdgeObjectIndexByIndexFrom(v,j);
            for ( int k = 0; k < graph.FromSize(v); ++k ) {
                if ( j == k ) continue;
                if ( graph.From(v)[k] == w ) continue;
                if ( pathsBetweenWithout( graph, v, w, edgei ) > 0 )
                    to_delete.push_back(j);
            }
        }
        UniqueSort(to_delete);
        graph.DeleteEdgesFrom(v,to_delete);
    }

    // ambiguity cleanup
    for ( int v = 0; v < graph.N(); ++v ) {
        if ( graph.FromSize(v) > 1 ) {
            vec<int> jvec, wvec;
            vec<int> j_to_delete;
            for ( int j = 0; j < graph.FromSize(v); ++j ) {
                jvec.push_back(j);
                wvec.push_back(graph.From(v)[j]);
            }
            SortSync(wvec,jvec);
            int nspecies = 1;
            for ( int i = 1; i < wvec.isize(); ++i )
                if ( wvec[i] != wvec[i-1] ) nspecies++;
            if ( wvec.size() > 2 && nspecies == 2 &&
                    wvec[0] != wvec[1] )
                j_to_delete.push_back(0);
            else if ( wvec.size() > 2 && nspecies == 2 &&
                    wvec[wvec.size()-1] != wvec[wvec.size()-2] )
                j_to_delete.push_back(jvec.size()-1);
            else if ( nspecies > 1 )
                j_to_delete = jvec;

            UniqueSort(j_to_delete);
            graph.DeleteEdgesFrom(v,j_to_delete);
        }
        if ( graph.ToSize(v) > 1 ) {
            vec<int> jvec, wvec;
            vec<int> j_to_delete;
            for ( int j = 0; j < graph.ToSize(v); ++j ) {
                jvec.push_back(j);
                wvec.push_back(graph.To(v)[j]);
            }
            SortSync(wvec,jvec);
            int nspecies = 1;
            for ( int i = 1; i < wvec.isize(); ++i )
                if ( wvec[i] != wvec[i-1] ) nspecies++;
            if ( wvec.size() > 2 && nspecies == 2 &&
                    wvec[0] != wvec[1] )
                j_to_delete.push_back(0);
            else if ( wvec.size() > 2 && nspecies == 2 &&
                    wvec[wvec.size()-1] != wvec[wvec.size()-2] )
                j_to_delete.push_back(wvec.size()-1);
            else if ( nspecies > 1 )
                j_to_delete = jvec;

            UniqueSort(j_to_delete);
            graph.DeleteEdgesTo(v,j_to_delete);

        }
    }


}



#if 0
void createLineInvolution( LineInfoVec const& lineInfo, vec<int>& lineInv )
{
    lineInv.clear();
    lineInv.resize( lineInfo.size() );
    for ( int line = 0; line < lineInfo.isize(); ++line ) {
        if ( lineInfo[line].mLineIdRC >= 0 ) {
            lineInv[line] = lineInfo[line].mLineIdRC;
        }
    }
}
#endif


void initializeLinesFromLoci( vec<LineMapReader::LineLocus> const& lineMap,
        int START_ID, int NUM_LINES, LineInfoVec& liv, vec<int>& lines )
{
    lines.clear();

    auto itr = lineMap.begin();

    if ( START_ID != -1 ) {
        while ( itr != lineMap.end() && itr->lineId != START_ID ) itr++;
        if ( itr == lineMap.end() ) {
            FatalErr("specified START_ID was not found");
        }
    }

    int num_lines = NUM_LINES;
    for ( /* see above */; itr != lineMap.end() && num_lines--; ++itr ) {
        auto mapLine = *itr;
        if ( mapLine.ref_gap || mapLine.len <= 6000 ) continue;
        int seed = mapLine.lineId;
        ForceAssertLt(seed, liv.isize());
        if ( liv[seed].mIsRC ) seed = liv[seed].mLineIdRC;
        lines.push_back(seed);
    }
}

void initializeScaffoldsFromLoci( vec<LineMapReader::LineLocus> const& lineMap,
        int START_ID, int NUM_LINES, int MINLINESIZE, vec<int> lineInv, ScaffoldPile& scaff )
{
    scaff.clear();

    auto itr = lineMap.begin();

    if ( START_ID != -1 ) {
        while ( itr != lineMap.end() && itr->lineId != START_ID ) itr++;
        if ( itr == lineMap.end() ) {
            FatalErr("specified START_ID was not found");
        }
    }

    int num_lines = NUM_LINES;
    for ( /* see above */; itr != lineMap.end() && num_lines--; ++itr ) {
        auto mapLine = *itr;
        if ( mapLine.ref_gap || mapLine.len <= MINLINESIZE ) continue;
        int seed = mapLine.lineId;
        ForceAssertLt(seed, lineInv.isize());
#warning "not canonicalizing for the moment"
//        if ( lineInv[seed] != -1 && lineInv[seed] < seed )
//            seed = lineInv[seed];
        vec<int> tmp = {seed };
        scaff.add(tmp);
    }
}


void linearComponents(LineGraph& graph, vec<vec<int>>& components)
{
    // This is simply for enumerating connected components of a guaranteed
    // linear (e.g. by construction) graph.  Things will break otherwise.
    vec<vec<int>> connected;
    graph.Components(connected);

    components.clear();
    for (int i = 0; i < connected.isize(); ++i) {
        for ( int j = 0; j < connected[i].isize(); ++j )
            if ( graph.Source(connected[i][j]) ) {
                vec<int> component;
                int v = connected[i][j];
                while (1) {
                    component.push_back(v);
                    if (graph.FromSize(v) == 0 ) break;
                    int from = graph.From(v)[0];
                    for ( int k = 1; k < graph.FromSize(v); ++k )
                        ForceAssertEq(from,graph.From(v)[k]);
                    v = from;
                }
                components.push_back(component);
                break;
            }
    }
}

};



int main( int argc, char* argv[] )
{
    RunTime( );

    BeginCommandArguments;
    CommandArgument_String_Doc(ADIR,"assembly dir.");
    CommandArgument_String_Doc(HICDIR,"hic.* files from CalcHiCLineHits");
    CommandArgument_String_Doc(LINEPAIRS,"Hi-C LINE-pairing file.");
    CommandArgument_String_Doc(SEPDIST,"hic.*.sepdist file");
    CommandArgument_Int(DEPTH);
    CommandArgument_Int_OrDefault_Doc(START_ID, -1, "starting line ID");
    CommandArgument_Int_OrDefault_Doc(NUM_LINES, -1, "starting line ID");
    CommandArgument_Bool_OrDefault(VERBOSE,True);
    CommandArgument_String_OrDefault_Doc(OUTHEAD,"", "output graph head" );
    CommandArgument_String_OrDefault_Doc(INHEAD,"", "input graph head (start mid-stream)" );
    CommandArgument_String_OrDefault_Doc(FINALHEAD,"", "final output graph head" );
    CommandArgument_Bool_OrDefault(FILTERSELF,True);
    CommandArgument_Int_OrDefault_Doc(MINLINESIZE,6000,"ignore lines smaller than this");
    CommandArgument_IntSet_OrDefault_Doc(TRACE_SEEDS,"","seeds to trace");
    CommandArgument_Double_OrDefault(VTHRESH,1e-5);
    CommandArgument_Int_OrDefault_Doc(MAXITER, std::numeric_limits<unsigned>::max(), "number of iterations to run");
    EndCommandArguments;

    String mapFile = ADIR + "/" + "a.lines.map";
    cout << Date() << ": reading line map from " << mapFile << endl;
    auto lineMap = LineMapReader(mapFile).LineLoci();
    cout << Date() << ": read " << lineMap.size() << " line alignments." << endl;

    HiCExperiment hic(
            ADIR + "/a.hbx",          // these are currently being loaded for no reason...
            ADIR + "/a.inv",
            HICDIR + "/" + EdgeIdToLineIdMap::gFileName,
            HICDIR + "/" + LineInfo::gFileName,
            ADIR + "/" + "a.lines",
            HICDIR + "/hic.hitrates",
            SEPDIST
            );

    // stash away line lengths
    // note that this could just be a vec... a relic right now that needs
    // to be done away with.
    unordered_map<int,int> lineLengths;
    for ( int lineId = 0; lineId < hic.mLineInfoVec.isize(); ++lineId )
        lineLengths[lineId] = hic.mLineInfoVec[lineId].mLineLen;


    cout << Date() << ": creating line involution" << endl;
    vec<int> lineInv;
    createLineInvolution(hic.mLineInfoVec, lineInv);

    cout << Date() << ": loading linePairs from " << LINEPAIRS << endl;
    LineLocPairVec linePairs;
    BinaryReader::readFile(LINEPAIRS, &linePairs);
    cout << Date() << ": " << linePairs.size() << " pairs loaded." << endl;

    /////// TEMP CHECK
    for ( auto const& llp : linePairs ) {
        if ( hic.mLineInfoVec[llp.ll1().getLineId()].mIsRC )
            cout << "L1 is RC" << endl;
        if ( hic.mLineInfoVec[llp.ll2().getLineId()].mIsRC )
            cout << "L2 is RC" << endl;
    }

    if ( FILTERSELF ) {
        linePairs.erase( std::remove_if( linePairs.begin(), linePairs.end(),
                [](LineLocPair const& llp) -> bool {
                return (llp.ll1().getLineId() == llp.ll2().getLineId());
        } ), linePairs.end() );
        linePairs.shrink_to_fit();
        cout << Date() << ": " << linePairs.size()  << " after filtering." << endl;
    }


    ScaffoldPile scaff(lineInv);
    initializeScaffoldsFromLoci( lineMap, START_ID, NUM_LINES, MINLINESIZE, lineInv, scaff );
//    for ( int i = 0; i < lineInv.isize(); ++i ) {
//        auto cline = canonicalLine(i, hic.mLineInfoVec);
//        if ( !scaff.has_key(cline) && hic.mLineInfoVec[i].mLoc.getChrId() == 1) {
//            scaff.add({cline});
//        }
//    }


    std::map<TaggedLineId, int> mapLineToVertex;
    int num_lines = NUM_LINES;

    cout << "initial scaffold: ";
    for ( auto itr = scaff.begin(); itr != scaff.end(); ++itr )
        cout << itr->first << ",";
    cout << endl;


    Bool changed = True;

    int iteration = 0;
    LineGraph graph;    // EdgeEvidence unused here
    while ( changed && iteration < MAXITER ) {
        graph.Clear();
        map<TaggedLineId,int> mapLineToVertex;

        changed = False;
        iteration++;

        set<int> seen;
        size_t count = 0;

        vec<int> seeds;
        for ( auto itr = scaff.begin(); itr != scaff.end(); ++itr ) {
            seeds.push_back(itr->first);
#ifdef GRAPH
            graph.AddVertex(TaggedLineId(itr->first, false));
            graph.AddVertex(TaggedLineId(itr->first, true));
#endif
        }
//        vec<int> seeds = scaff.starts();

        for ( auto seed : seeds ) {
            count++;

            if ( seen.find(seed) != seen.end() ) continue;

            PRINT3(iteration, count, seed);

            Bool verbose = VERBOSE && Member(TRACE_SEEDS,seed);

            auto agg = hiCGrowNeighborhoodScaff(seed, scaff, linePairs,
                    lineInv, lineLengths,
                    hic, DEPTH, MINLINESIZE, iteration, verbose, OUTHEAD, VTHRESH);

            cout << "seed=" << seed << ", agg=" << agg[0] << "," << agg[1] << "," << agg[2] << endl;

            // if somehow left and right are the same scaffold, we bail
            // on this endeavor
            if ( agg[0].getId() != TaggedLineId::maxVal &&
                    agg[2].getId() != TaggedLineId::maxVal &&
                    scaff.start(agg[0].getId()) == scaff.start(agg[2].getId()))
                continue;

            vec<int> to_merge;
            vec<int> sizes;
            if ( agg[0].getId() != TaggedLineId::maxVal ) {
                // there was a left, so try to grow it
                if ( verbose ) cout << "============== aggl for seed "
                        << seed << " =================" << endl;
                auto aggl = hiCGrowNeighborhoodScaff(agg[0].getId(), scaff,
                        linePairs, lineInv, lineLengths, hic, DEPTH, MINLINESIZE, iteration, verbose);
                // okay, so make sure that the id wasn't canonicalized --
                // it already should be in the scaffold
                ForceAssertEq(aggl[1].getId(), agg[0].getId());
                // and the results should be relative to its forward orientation
                ForceAssertEq(aggl[1].isTag(), false);
                // but if our 'agg' result showed it flipped, then flip this
                // result
                if ( agg[0].isTag() ) flipTaggedArr(aggl);
                // okay, now we should have precise equality on the seed in aggl
                ForceAssert( aggl[1] == agg[0] );
                cout << "seed=" << seed << ", aggl=" << aggl[0] << ","
                        << aggl[1] << "," << aggl[2] << endl;
                if ( aggl[2] == agg[1] ) {
#ifdef GRAPH
                    addToLineGraph(agg[0], agg[1], graph, mapLineToVertex);
                    changed=True;
#else
                    // okay, accept the left side
                    // rc the scaffold, if necessary
                    if ( agg[0].isTag() ) scaff.rc(agg[0].getId());
                    to_merge.push_back( agg[0].getId() );
                    sizes.push_back( scaff.scaffold(to_merge.back()).size() );
#endif
                }
            }
            to_merge.push_back( agg[1].getId() );
            sizes.push_back( scaff.scaffold(to_merge.back()).size() );
            if ( agg[2].getId() != TaggedLineId::maxVal ) {
                if ( verbose ) cout << "============== aggr for seed "
                        << seed << " =================" << endl;
                auto aggr = hiCGrowNeighborhoodScaff(agg[2].getId(), scaff,
                        linePairs, lineInv, lineLengths, hic, DEPTH, MINLINESIZE, iteration, verbose);
                ForceAssertEq(aggr[1].getId(), agg[2].getId());
                ForceAssertEq(aggr[1].isTag(), false);

                if ( agg[2].isTag() ) flipTaggedArr(aggr);
                // okay, now we should have precise equality on the seed in aggr
                ForceAssert(aggr[1] == agg[2]);
                cout << "seed=" << seed << ", aggr=" << aggr[0] << ","
                        << aggr[1] << "," << aggr[2] << endl;
                if ( aggr[0] == agg[1] ) {
#ifdef GRAPH
                    addToLineGraph(agg[1], agg[2], graph, mapLineToVertex);
                    changed=True;
#else
                    // okay, accept the right side
                    if ( agg[2].isTag() ) scaff.rc(agg[2].getId());
                    to_merge.push_back( agg[2].getId() );
                    sizes.push_back( scaff.scaffold(to_merge.back()).size() );
#endif
                }
            }

            cout << "to_merge: " << printSeq(to_merge) << ", sizes: " << printSeq(sizes) << endl;
            if ( 1 ) {
                for ( auto m : to_merge ) {
                    cout << m  << " inverse is " << lineInv[m] << endl;
                    for ( auto s : scaff.scaffold(m) ) {
                        seen.insert(s);
                        seen.insert(lineInv[s]);
                    }
                }
            }
//            if (  to_merge.size() > 2 || ( to_merge.size() > 1 && !Member(sizes,1) ) ) {
            if (  to_merge.size() > 1 ) {
                changed=True;
                cout << "merging: " << printSeq(to_merge) << endl;
                for ( auto line : to_merge ) cout << line << ": "
                        << printSeq(scaff.scaffold(line)) << endl;
                scaff.merge(to_merge);
                cout << "merged: " << printSeq(scaff(to_merge[0])) << endl;
            }
        }


#ifdef GRAPH
        dumpLineGraph(graph, OUTHEAD+".dirty.dot", -1, False, False);
        cleanAmbiguities(graph, hic, VERBOSE );
        dumpLineGraph(graph, OUTHEAD+".clean.dot", -1, False, False);
        // see if the start of each scaffold is connected to something else and
        // merge, if necessary
        vec<int> starts = scaff.starts();
        for ( int i = 0; i < graph.N(); ++i )
            if ( Member(starts, graph.Vert(i) ) && graph.ToSize(i) == 1 ) {
                auto toId = graph.Vert(graph.To(i)[0]);
            }
#endif


    }

#ifdef GRAPH
    BinaryWriter::writeFile(OUTHEAD + ".graph", graph);

//    cleanAmbiguities(graph, hic, VERBOSE );

    vec<vec<int>> components;
    linearComponents(graph,components);

    // turn line graph into scaffolds and let's see how we did
    std::map<TaggedLineId,pair<int,int>> mapLineToScaffold;
    for ( int i = 0; i < components.isize(); ++i ) {
        if ( components[i].size() == 1 ) continue;
        for ( int j = 0; j < components[i].isize(); ++j ) {
            auto v = components[i][j];
            mapLineToScaffold[graph.Vert(v)] = std::make_pair(i,j);
        }
    }
#else
    std::map<TaggedLineId,pair<int,int>> mapLineToScaffold;
    vec<vec<int>> components = scaff.scaffolds();
    for ( int i = 0; i < components.isize(); ++i ) {
        if ( components[i].size() == 1 ) continue;
        for ( int j = 0; j < components[i].isize(); ++j ) {
            int line = components[i][j];
            int cline = canonicalLine(line, hic.mLineInfoVec);
            bool rc = (line != cline);
            mapLineToScaffold[ TaggedLineId(cline, rc) ] = make_pair(i,j);
        }
    }
#endif

    auto itr = lineMap.begin();
    if ( START_ID != -1 ) {
        while ( itr != lineMap.end() && itr->lineId != START_ID ) itr++;
        if ( itr == lineMap.end() ) {
            FatalErr("specified START_ID was not found");
        }
    }

    int last_scaf = -1, last_line = -1;
    bool last_rc = false;
    bool first = true;
    for ( int num_lines = NUM_LINES; itr != lineMap.end() && num_lines--; ++itr ) {
        auto mapLine = *itr;
//        if (mapLine.ref_gap || mapLine.len <= MINLINESIZE ) continue;
        if ( mapLine.ref_gap ) continue;
        String mark=" ";
        printCols(3, mark);
        ostringstream ourStr;

        if ( mapLine.len > MINLINESIZE ) {
            TaggedLineId forw( mapLine.lineId, false );
            TaggedLineId rev( mapLine.lineId, true );

            bool found_forw = mapLineToScaffold.find( forw ) != mapLineToScaffold.end();
            bool found_rev = mapLineToScaffold.find( rev ) != mapLineToScaffold.end();

            int this_scaf=-1, this_line=-1;
            bool this_rc=false;

            if ( !found_forw && !found_rev ) mark = "*";

            if ( found_forw ) {
                this_scaf = mapLineToScaffold[forw].first;
                this_line = mapLineToScaffold[forw].second;
                this_rc = false;
                ourStr << "+" << this_scaf << "." << this_line;
            }
            if ( found_forw && found_rev ) {
                ourStr << ",";
                mark = "*";
            }
            if ( found_rev ) {
                this_scaf = mapLineToScaffold[rev].first;
                this_line = mapLineToScaffold[rev].second;
                this_rc = true;
                ourStr << "-" << this_scaf << "." << this_line;
            }

            if (first && this_scaf != -1)
                first = false;
            else if ( this_scaf != -1 ) {
                if ( this_scaf != last_scaf ) {    // did we end properly?
                    if ( last_rc == true && last_line != 0 ) mark="*";
                    else if ( last_rc == false && last_line
                            != components[last_scaf].isize() - 1 )
                        mark = "*";
                    else if ( this_rc == false && this_line != 0 )
                        mark = "*";
                    else if ( this_rc == true && this_line != components[this_scaf].isize()-1 )
                        mark = "*";
                }
                else if ( this_rc != last_rc ) mark="*";
                else if (!this_rc && this_line != last_line+1 ) mark="*";
                else if (this_rc && this_line != last_line-1 ) mark="*";
            }

            if ( this_scaf != -1 ) {
                last_scaf = this_scaf;
                last_rc = this_rc;
                last_line = this_line;
            }
        }

        printCols(3,mark);
        printCols(20, ourStr.str());

        printCols(20, mapLine.chr + ":" + ToString(mapLine.start) +
                "-" + ToString(mapLine.stop));
        printCols(14, " [" + ToString(mapLine.lineId) + "]: " );

        ostringstream lrs;
        printCols(20,lrs.str());
        printCols(15," len=" + ToString(mapLine.len) + " cov=" + ToString(mapLine.cov));
        cout << endl;
    }


//    hiCSuperScaffold(graph, hic, DEPTH, VERBOSE, OUTHEAD);

    if (FINALHEAD != "") {
#ifdef GRAPH
        dumpLineGraph(graph, FINALHEAD+".final.dot", -1, False, False);
#endif
        Ofstream(idx, FINALHEAD+".final.txt");
        for ( size_t i = 0; i < components.size(); ++i ){
            idx << "[scaff " <<  i << "]: ";
            for ( auto const l : components[i] )
                idx << l << " ";
            idx << endl;
        }
    }

#if 0
    for ( ; itr != lineMap.end(); ++itr ) {

        if ( NUM_LINES > 0 && --NUM_LINES == 0 ) break;
    }
#endif


    return 0;
}



