///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * ReadPairImpliedAdjacencies.cc
 *
 *  Created on: Dec 4, 2013
 *      Author: neilw
 *
 * This was created specifically for a data-type where read pairs imply
 * potentially long distance relationships.
 *
 */
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
#include <omp.h>
#include <iomanip>
#include "MainTools.h"
#include "PairsManager.h"
#include "paths/HyperBasevector.h"
#include "paths/long/LocTraceEdges.h"
#include "paths/long/EvalByReads.h"
#include "paths/UnibaseUtils.h"
#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"
#include <queue>

struct Adjacency {
    size_t edge2;
    size_t read1;
    size_t read2;
};
typedef SerfVec<Adjacency> AdjacencyVec;
typedef MasterVec<AdjacencyVec> VecAdjacencyVec;
// feudal vector template instantiation at end of file... will
// break into header/implementation pair later... next line will
// end up in the header.
//extern template class OuterVec<AdjacencyVec>;

bool PerfectReadPlace( read_place const& place, HyperBasevector const& hbv,
	basevector const& bases, bool debug = false )
{
    size_t edge0 = place.E(0);
    basevector bedge0 = hbv.EdgeObject(edge0);
    if (!place.Fw()) bedge0.ReverseComplement();
    size_t matchlen = min( bedge0.size() - place.P(), bases.size() );
    basevector bv1, bv2;
    bv1.SetToSubOf(bedge0,place.P(),matchlen);
    bv2.SetToSubOf(bases, 0, matchlen);

    if ( debug ) {
	cout << "DEBUG READ" << endl;
	cout << "EDGE:" << endl;
	bv1.Print(cout);
	cout << "READ:" << endl;
	bv2.Print(cout);
	cout << "MATCH: " << (bv1 == bv2) << endl;
    }

    return (bv1 == bv2);
}


VecAdjacencyVec AlignPairs2( const HyperBasevector& hbv, const vecbasevector& bases,
	const vecqualvector& quals, PairsManager const& pairs,
        vecString* namesp = nullptr,
	const bool require_perfect = false, const bool debug = false )
{
    int L = 12;
    double DELTA_QSUM = 0.0;

    vec<vec<read_place>> places_log = LocTraceEdges( hbv, bases, quals, DELTA_QSUM, L, false );

    VecAdjacencyVec adj( hbv.EdgeObjectCount() );

    cout << "TABLE:" << endl;

    for ( size_t i = 0; i < pairs.nPairs(); ++i ) {
	size_t id1 = pairs.ID1(i);
	size_t id2 = pairs.ID2(i);
	auto const& places1 = places_log[id1];
	auto const& places2 = places_log[id2];
	if ( places1.size() == 1 && places2.size() == 1 ) {
	    // note that we're just taking the first alignment if there are multiple
	    // edges to which the read aligns
	    if ( !require_perfect || ( places1[0].Qsum() == 0 && places2[0].Qsum() == 0 ) ) {
		size_t edge1 = places1[0].E(0);
		size_t edge2 = places2[0].E(0);
		adj[edge1].push_back( Adjacency{ edge2, id1, id2 } );
		if ( edge1 != edge2 ) adj[edge2].push_back( Adjacency{ edge1, id1, id2 } );
                //NIW
                if ( namesp ) {
                    vecString const& names = *namesp;
                    ForceAssertEq(names[id1],names[id2]);
                    cout << names[id1] << endl;
                }
	    }
	}
    }

    cout << "END OF TABLE" << endl;

    for_each( adj.begin(), adj.end(), [](AdjacencyVec& v) {
	std::sort( v.begin(), v.end(), []( Adjacency const& v1, Adjacency const& v2 ) {
	    return v1.edge2 < v2.edge2; });
    } );

    if ( 0 && debug ) {
        size_t pair_placements = 0;
	for ( size_t i = 0; i < pairs.nPairs(); ++i ) {
	    size_t id1 = pairs.ID1(i);
	    size_t id2 = pairs.ID2(i);
	    cout << "pair " << i << ": reads (" << id1 << "," << id2 << ") ----------------------------" << endl;
	    cout << "\tleft:" << endl;
	    cout << "\t\t";
	    auto const& places0=places_log[id1];
	    std::copy( places0.begin(), places0.end(),
		    std::ostream_iterator<read_place>(cout, "\n\t\t"));
	    cout << endl;

            cout << "READ ID #" << id1 << endl;
            bases[id1].Print(cout);
	    if ( places0.size() ) {
		size_t edge0 = places0[0].E(0);
		basevector bedge0 = hbv.EdgeObject(edge0);
		if (!places0[0].Fw() ) bedge0.ReverseComplement();
		cout << "EDGE #" << edge0 << endl;
		bedge0.PrintBases(cout, places0[0].P(), bases[id1].size() );
	    } else {
                cout << "NO PLACEMENTS" << endl;
            }

	    cout << endl << "\tright:" << endl;
	    cout << "\t\t";
	    auto const& places1=places_log[id2];
	    std::copy( places1.begin(), places1.end(),
		    std::ostream_iterator<read_place>(cout, "\n\t\t"));
	    cout << endl;

            cout << "READ ID #" << id2 << endl;
            bases[id2].Print(cout);
	    if ( places1.size() ) {
		size_t edge0 = places1[0].E(0);
		basevector bedge0 = hbv.EdgeObject(edge0);
		if (!places1[0].Fw()) bedge0.ReverseComplement();
		cout << "EDGE #" << edge0 << endl;
		bedge0.PrintBases(cout,places1[0].P(), bases[id2].size());
	    } else {
                cout << "NO PLACEMENTS" << endl;
            }

            if ( places0.size() && places1.size() )
                pair_placements++;
	}
        cout << Date() << ": found placements for " << pair_placements
            << " of " << pairs.nPairs() << " pairs." << endl;
    }


    return adj;
}

vec<int> find_bubble_edges( HyperBasevector const& hbv, size_t v1, size_t v2 )
{
    vec<int> bubble_edges;
    // find all of the bubble edges, in order
    size_t currv = v1;
    while ( currv != v2 ) {
        if ( hbv.FromSize(currv) == 2 )
            bubble_edges.push_back(
                hbv.FromEdgeObj(currv)[0],
                hbv.FromEdgeObj(currv)[1] );
        ForceAssertGe( hbv.FromSize(currv), 1 );
        currv = hbv.From(currv)[0];
    }
    return bubble_edges;
}

vec<int> count_evidence_in_set(VecAdjacencyVec const& adj,
            vec<int> const& bubble_edges )
{
    ForceAssertEq( bubble_edges.size() % 2 , 0U );
    vec<int> counts( bubble_edges.size(), 0);
    for ( size_t i = 0; i < bubble_edges.size(); ++i ) {
        size_t edgeno = bubble_edges[i];
        size_t other_edgeno = ( i % 2 == 0 ) ? bubble_edges[i+1] : bubble_edges[i-1];
        for ( auto const& a : adj[edgeno] ) {
            if ( Member( bubble_edges, (int)a.edge2 ) && a.edge2 != edgeno &&
                    a.edge2 != other_edgeno )   // pairs to an edge of ANOTHER bubble
                counts[i]++;
        }
    }
    return counts;
}


void part_bubble_runs_on_evidence( HyperBasevector const& hbv,
        VecAdjacencyVec const&adj, std::vector< std::pair<size_t, size_t> >& runs,
        vec<int>& parts,
        bool debug = false )
{
    std::vector<std::pair<size_t,size_t> > new_runs;

    for ( auto const run : runs ) {
        auto bubble_edges = find_bubble_edges( hbv, run.first, run.second );

        // compute scoring matrix
        size_t nedges = hbv.EdgeObjectCount();
        vec<vec<int>> adj_adj(nedges);
        for ( size_t i = 0; i < nedges; ++i ) {
            if ( !Member( bubble_edges, (int)i ) ) continue;
            adj_adj[i].resize(nedges,0);
            for ( auto const& a : adj[i] )
                if ( Member( bubble_edges, (int)a.edge2 ) )
                    adj_adj[i][a.edge2]+=1;
        }

        int this_part = 0, new_part = 0;
        parts.resize( hbv.EdgeObjectCount(), 0 );
        vec<bool> seen( parts.size(), 0 );
        std::map<int,int> comp_edge;
        for ( size_t i = 0; i < bubble_edges.size(); i+=2 ) {
            int e1 = bubble_edges[i];
            int e2 = bubble_edges[i+1];
            comp_edge[e1] = e2;
            comp_edge[e2] = e1;
        }
        vec<int> comp_parts( bubble_edges.size() + 2 );
        for ( size_t i = 2; i < bubble_edges.size()+2; i+=2 ) {
            comp_parts[i] = i+1;
            comp_parts[i+1]=i;
        }



        auto bubble_edge_copy( bubble_edges );

	while (bubble_edges.size()) {
	    int edge = bubble_edges.back();
	    bubble_edges.pop_back();
	    if (seen[edge])
		continue;
	    else
		seen[edge] = true;
	    if (seen[comp_edge[edge]])
		this_part = comp_parts[parts[comp_edge[edge]]];
	    else
		this_part = (new_part += 2);
	    seen[comp_edge[edge]] = true;
	    cout << "unvisited bubble edge " << edge << " gets " << this_part
		    << endl;
	    cout << "complementary edge " << comp_edge[edge] << " gets "
		    << comp_parts[this_part] << endl;
	    parts[edge] = this_part;
	    parts[comp_edge[edge]] = comp_parts[this_part];
	    cout << "THIS LABEL " << this_part << " based on edge " << edge
		    << endl;
	    vec<int> cur_part;
	    for (size_t i = 0; i < adj_adj[edge].size(); ++i)
		if (adj_adj[edge][i] > 0 && !seen[i])
		    cur_part.push_back(i);
	    while (cur_part.size()) {
		int this_edge = cur_part.back();
		cur_part.pop_back();
		if (seen[this_edge])
		    continue;
		else {
		    seen[this_edge] = true;
		    seen[comp_edge[this_edge]] = true;
		}
		PRINT2(edge, this_edge);
		std::copy(cur_part.begin(), cur_part.end(),
			std::ostream_iterator<int>(cout, " "));
		cout << "reached via " << edge << " bubble edge " << this_edge
			<< " gets " << this_part << endl;
		cout << "complementary edge " << comp_edge[this_edge]
			<< " gets " << comp_parts[this_part] << endl;
		parts[this_edge] = this_part;
		parts[comp_edge[this_edge]] = comp_parts[this_part];
		for (size_t i = 0; i < adj_adj[this_edge].size(); ++i)
		    if (adj_adj[this_edge][i] > 0 && !seen[i])
			cur_part.push_back(i);
	    }
	}

	cout << "GRAPH LABELING for run from v1=" << run.first << " to v2="
		<< run.second << endl;
	for (auto const& edge : bubble_edge_copy) {
	    cout << edge << " (part " << parts[edge] << ")" << endl;
	}

        size_t currv = run.first;
        size_t startv = run.first;
        size_t endv = run.first;
        size_t bubble_n = 0;
        vec<int> last;
        while ( currv != run.second ) {
            cout << "work: ";
            PRINT3( currv, startv, endv );
            if ( hbv.FromSize(currv) == 2 ) {
                int edge0 = hbv.FromEdgeObj(currv)[0];
                int edge1 = hbv.FromEdgeObj(currv)[1];
                if ( last.size() == 0 ) {
                    last.push_back( parts[edge0], parts[edge1] );
                }
                else if (!Member(last, parts[edge0]) ||
                	 !Member(last, parts[edge1])) {
                    PRINT7(currv, edge0, edge1, last[0], last[1], parts[edge0], parts[edge1] );
                    if ( endv != startv ) new_runs.push_back(std::make_pair( startv, endv) );
                    startv = currv;
                    last[0] = parts[edge0];
                    last[1] = parts[edge1];
                }
                bubble_n++;
                endv = hbv.From(currv)[0];
            }
            currv = hbv.From(currv)[0];
        }
        if ( endv != startv )
            new_runs.push_back(std::make_pair( startv, endv) );
    }
    swap(runs, new_runs);
}

void split_bubble_runs_on_evidence( HyperBasevector const& hbv,
        VecAdjacencyVec const&adj, std::vector< std::pair<size_t, size_t> >& runs,
        bool debug = false )
{
    std::vector< std::pair<size_t, size_t> > new_runs;

    for ( auto const run : runs ) {
        auto bubble_edges = find_bubble_edges( hbv, run.first, run.second );
        auto bubble_counts = count_evidence_in_set( adj, bubble_edges );

        std::cout << "BUBBLE_COUNTS=";
        std::copy( bubble_counts.begin(), bubble_counts.end(),
                std::ostream_iterator<int>( cout, " " )) ;
        std::cout << endl;

        size_t currv = run.first;
        size_t startv = run.first;
        size_t endv = run.second;
        size_t bubble_n = 0;
        while ( currv != run.second ) {
            if ( hbv.FromSize(currv) == 2 ) {
                PRINT5(startv, currv, bubble_n, bubble_counts[2*bubble_n], bubble_counts[2*bubble_n+1]);
                if ( bubble_counts[2*bubble_n] + bubble_counts[2*bubble_n+1] == 0 ) {
                    if ( endv != startv )
                        new_runs.push_back(std::make_pair( startv, endv) );
                    startv = hbv.From(currv)[0];
                }
                bubble_n++;
                endv = hbv.From(currv)[0];
            }
            currv = hbv.From(currv)[0];
        }
        if ( endv != startv )
            new_runs.push_back(std::make_pair( startv, endv) );
    }

    cout << "NEW_RUNS:" << endl;
    for ( auto const& run : new_runs ) {
        cout << "(" << run.first << "," << run.second << ") ";
    }
    cout << endl;

    std::swap( new_runs, runs );
}

std::vector< std::pair<size_t, size_t> >
find_bubble_runs( HyperBasevector const& hbv, bool debug = false )
{
    debug=true;
    std::vector<std::pair<size_t,size_t> > runs;

    std::vector<bool> visited( hbv.N() );
    vec<vec<int> > comps;
    hbv.Components( comps );

    for ( size_t ci = 0; ci < comps.size(); ++ci ) {
	if (debug) cout << "component #" << ci << endl;
	vec<int> sources;
	hbv.SubgraphSources( comps[ci], sources );
	if (debug) cout << "sources=" << sources << endl;

	for ( auto source : sources ) {
	    std::vector<size_t> to_visit;
	    bool in_bubble_run = false;

	    to_visit.push_back(source);

	    while ( to_visit.size() ) {

		size_t currv = to_visit.back();
		to_visit.pop_back();

		if ( visited[currv] )
		    in_bubble_run = false;	// avoid crazy bubble cycles
		else {
		    visited[currv] = true;

		    // define conditions for starting and continuing a bubble run
		    bool bubble_start = hbv.FromSize(currv) == 2 && hbv.From(currv)[0] == hbv.From(currv)[1] && hbv.ToSize(hbv.From(currv)[0]) == 2;
		    bool bubble_cont = bubble_start || ( hbv.FromSize(currv) == 1 && hbv.ToSize(hbv.From(currv)[0]) == 1 );

		    // push back next vertices to visit
		    to_visit.insert(to_visit.end(), hbv.From(currv).begin(), hbv.From(currv).end() );	// might double up; doesn't matter

		    // always mark the current vertex as a potential terminator for the
		    // bubble run in case it hits a cycle and gets cut off next iteration
		    if ( in_bubble_run ) runs.back().second=currv;

		    if ( in_bubble_run && !bubble_cont ) {
			in_bubble_run = false;
		    } else if (!in_bubble_run && bubble_start ) {
			in_bubble_run = true;
			runs.push_back( std::make_pair( currv, 0 ) );
		    }
		}

	    }

	}
    }

    // don't end a run on a single edge -- easier to do here
    for ( auto& run : runs ) {
	if ( hbv.ToSize(run.second) == 1 )
	    run.second = hbv.To(run.second)[0];
    }


    return runs;
}


vec<vec<int>> explode_bubble_edges( vec<int> const& bubble_edges )  {
    ForceAssertEq( bubble_edges.size() % 2, 0U);
    size_t nbubbles = bubble_edges.size() / 2;
    size_t new_paths = 1 << nbubbles;
    vec<vec<int>> kapow(new_paths);
    for ( size_t i = 0; i < new_paths; ++i ) {
	for ( size_t j = 0, temp=i; j < nbubbles; ++j, temp>>=1 )
	    kapow[i].push_back( bubble_edges[2*j+(temp&0x1U)] );
    }

    return kapow;
}

void edit_bubble_runs_based_on_parts( HyperBasevector& hbv,
	std::vector<std::pair<size_t,size_t> > const& runs,
	vec<int>& parts, bool debug = false)
{
    size_t orig_bubbles = 0;
    size_t new_bubbles = 0;

    auto const&  edges = hbv.Edges();
    vec<int> to_delete;
    for ( auto const& run : runs ) {
	orig_bubbles += find_bubble_edges(hbv, run.first, run.second).size() / 2;
	vec<int> last;
	basevector new_e0, new_e1;
	size_t currv = run.first;
	while ( currv != run.second ) {
	    if ( hbv.FromSize(currv) == 1 )  {
		cout << "side 1" << endl;
		size_t from0 = hbv.FromEdgeObj(currv)[0];
		if ( new_e0.size() ) new_e0 = TrimCat( hbv.K(), new_e0, edges[ from0 ] );
		else new_e0 = edges[from0];
		if ( new_e1.size() ) new_e1 = TrimCat( hbv.K(), new_e1, edges[ from0 ] );
		else new_e1 = edges[from0];
		to_delete.push_back( from0 );
	    } else if ( hbv.FromSize(currv) == 2 ) {
		cout << "side 2 " << endl;
		int from0 = hbv.FromEdgeObj(currv)[0];
		int from1 = hbv.FromEdgeObj(currv)[1];
		PRINT4(new_e0.size(), new_e1.size(), edges[from0].size(), edges[from1].size());
		// unnecessary parenthesis sponsored by Eclipse Indigo
		if ( last.size() == 0 || ( Member(last, parts[from0]) && Member(last, parts[from1]) ) )  {
		    if ( last.size() == 0 ) last.push_back( parts[from0], parts[from1] );
		    int edge0 = ( parts[from1] == last.back() ) ? from0 : from1;
		    int edge1 = ( parts[from1] == last.back() ) ? from1 : from0;
		    if ( new_e0.size() ) new_e0 = TrimCat( hbv.K(), new_e0, edges[ edge0 ] );
		    else new_e0 = edges[ edge0 ];
		    if ( new_e1.size() ) new_e1 = TrimCat( hbv.K(), new_e1, edges[ edge1 ] );
		    else new_e1 = edges[ edge1 ];
		    to_delete.push_back( from0 );
		    to_delete.push_back( from1 );
		} else {
		    FatalErr("BUG: lost synchromization with bubble edge list");
		}
	    }
	    currv = hbv.From(currv)[0];
	}
	hbv.AddEdge(run.first, run.second, new_e0);
	hbv.AddEdge(run.first, run.second, new_e1);
	hbv.TestValid();
	new_bubbles += 1;
    } // for runs...
    hbv.DeleteEdges(to_delete);
    hbv.RemoveDeadEdgeObjects();
    hbv.RemoveEdgelessVertices();

    cout << "RESULT: " << orig_bubbles << " to " << new_bubbles << endl;
}

void edit_bubble_runs( HyperBasevector& hbv,
	std::vector<std::pair<size_t,size_t> > const& runs,
	VecAdjacencyVec const& adj, bool debug = false)
{
    debug = true;

    // compute scoring matrix
    size_t nedges = hbv.EdgeObjectCount();
    vec<vec<int>> adj_adj(nedges);
    for ( size_t i = 0; i < nedges; ++i ) {
	adj_adj[i].resize(nedges,0);
	for ( auto const& a : adj[i] )
	    adj_adj[i][a.edge2]+=1;
    }

    // print scoring matrix
    cout << "SCORING MATRIX:" << endl;
    for ( size_t i = 0; i < nedges; ++i ) {
        cout << i << ": ";
        for ( size_t j = 0; j < nedges; ++j )
            if ( adj_adj[i][j] > 0 )
                cout << j << "(" << adj_adj[i][j] << ") ";
        cout << endl;
    }
    cout << "END OF SCORING MATRIX" << endl;


    // debugging: histogram coverage of edges
    vec<size_t> hist(nedges,0U);
    PRINT2( adj.size(), nedges );
    for ( size_t i = 0; i < nedges; ++i )
        for ( auto const& a : adj[i] )
            hist[ a.edge2 ]++;
    cout << "HISTOGRAM:" << endl;
    for ( size_t i = 0; i < nedges; ++i )
        cout << "hist[" << i << "]=" << hist[i] << endl;


    vec<int> to_delete;
    for ( auto const& run : runs ) {
	// okay, at this point, this is just an expensive sanity check
	vec<vec< int>> paths;
	const int maxpaths = 1024;
	if ( !hbv.AllPaths(run.first, run.second, paths, 1024) ) {
	    cout << "failed to expand AllPaths between " << run.first <<
		    "," << run.second << endl;
	    continue;
	}
	ForceAssertEq(paths.size(), 1U);

	// find all of the bubble edges, in order
	vec<int> bubble_edges;
	size_t currv = run.first;
	while ( currv != run.second ) {
	    if ( hbv.FromSize(currv) == 2 )
		bubble_edges.push_back(
		    hbv.FromEdgeObj(currv)[0],
		    hbv.FromEdgeObj(currv)[1] );
	    ForceAssertGe( hbv.FromSize(currv), 1 );
	    currv = hbv.From(currv)[0];
	}

	if ( 1|| debug ) {
	    cout << "edge run from vertex " << run.first << " to " << run.second << endl;
	    cout << "we found " << bubble_edges.size() << " bubble edges:" << endl;
	    std::copy( bubble_edges.begin(), bubble_edges.end(),
		    std::ostream_iterator<int>( cout, ",") );
	    cout << endl;

	}

        // print scoring matrix
        cout << "SCORING MATRIX FOR BUBBLE EDGES ONLY:" << endl;
        int iedges = static_cast<int>(nedges);
        for ( int i = 0; i < iedges; ++i ) {
            if (!Member( bubble_edges, i) ) continue;
            cout << i << ": ";
            for ( int j = 0; j < iedges; ++j ) {
                if (!Member( bubble_edges, j) ) continue;
                if ( adj_adj[i][j] > 0 )
                    cout << j << "(" << adj_adj[i][j] << ") ";
            }
            cout << endl;
        }
        cout << "END OF SCORING MATRIX FOR BUBBLE EDGES ONLY" << endl;


	vec< vec< int>> bubble_paths = explode_bubble_edges( bubble_edges );

	// make sparse adj_adj just for this path
	// oy -- ugly, but let's see what the data says before
	// we worry
	vec< std::map< size_t, size_t > > adj_map( nedges );
	for ( auto const& edge1 : bubble_edges )
	    for ( auto const& edge2 : bubble_edges )
		adj_map[edge1][edge2] = adj_adj[edge1][edge2];


	cout << "BUBBLE PATHS (" << bubble_paths.size() << "):" << endl;

	if ( debug ) {
            for ( size_t i = 0; i < bubble_paths.size(); ++i ) {
                auto const& path = bubble_paths[i];
                cout << i << ": ";
		std::copy( path.begin(), path.end(),
			std::ostream_iterator<int>( cout, "->") );
		cout << endl;
	    }
	}

	vec<int> score( bubble_paths.size(), 0);
	for ( size_t i = 0; i < bubble_paths.size(); ++i ) {
	    vec<int> const& path = bubble_paths[i];
	    for ( size_t j = 0; j < path.size(); ++j ) {
		auto const& adj_line_map = adj_map[path[j]];
		size_t score_delta = std::accumulate(
			adj_line_map.begin(),
			adj_line_map.end(), 0U,
			[](size_t p1, std::pair<size_t,size_t> p2){ return p1+p2.second;});

#if 0
                for ( auto const& m : adj_line_map ) {
                    std::cout << "(" << m.first << "," << m.second << ")" << endl;
                }
                std::cout << "=" << score_delta << endl;
#endif

                if ( debug ) cout << i << ": " << score_delta;

		for ( size_t k = 0; k < path.size(); ++k )
		    score_delta -= adj_line_map.at(path[k]);

                if ( debug ) cout << " (" << score_delta << ")" << endl;


                PRINT3( i, j, score_delta );
		score[i] += score_delta;
	    }
	}

	cout << "SCORES (sorted):" << endl;
        vec<int> scorei( score.size(), vec<int>::IDENTITY );
	SortSync(score,scorei);
        for ( size_t i = 0; i < score.size(); ++i )
            cout << score[i] << " index " << scorei[i] << endl;


        if ( score.size() > 1 && score[0] == 0 && score[1] == 0 && score[2] != 0 ) {
            vec<int> rpaths0, rpaths1;
            auto const& paths0 = bubble_paths[scorei[0]];
            auto const& paths1 = bubble_paths[scorei[1]];
            rpaths0.insert( rpaths0.end(), paths0.rbegin(), paths0.rend() );
            rpaths1.insert( rpaths1.end(), paths1.rbegin(), paths1.rend() );
            cout << "PATHS0: ";
            std::copy( paths0.begin(), paths0.end(), std::ostream_iterator<int>(cout, " ") );
            cout << endl;
            cout << "RPATHS0: ";
            std::copy( rpaths0.begin(), rpaths0.end(), std::ostream_iterator<int>(cout, " ") );
            cout << endl;
            cout << "RPATHS1: ";
            std::copy( rpaths1.begin(), rpaths1.end(), std::ostream_iterator<int>(cout, " ") );
            cout << endl;

            hbv.TestValid();

            basevector new_e0, new_e1;
            auto const&  edges = hbv.Edges();
            size_t currv = run.first;
            while ( currv != run.second ) {
        	if ( hbv.FromSize(currv) == 1 )  {
        	    cout << "side 1" << endl;
        	    size_t from0 = hbv.FromEdgeObj(currv)[0];
        	    if ( new_e0.size() ) new_e0 = TrimCat( hbv.K(), new_e0, edges[ from0 ] );
        	    else new_e0 = edges[from0];
        	    if ( new_e1.size() ) new_e1 = TrimCat( hbv.K(), new_e1, edges[ from0 ] );
        	    else new_e1 = edges[from0];
        	    to_delete.push_back( from0 );
        	} else if ( hbv.FromSize(currv) == 2 ) {
        	    cout << "side 2 " << endl;
        	    int from0 = hbv.FromEdgeObj(currv)[0];
        	    int from1 = hbv.FromEdgeObj(currv)[1];
        	    PRINT4(new_e0.size(), new_e1.size(), edges[from0].size(), edges[from1].size());
        	    // unnecessary parenthesis sponsored by Eclipse Indigo
        	    if ( ( from0 == rpaths0.back() && from1 == rpaths1.back() ) ||
        		    ( from0 == rpaths1.back() && from1 == rpaths0.back() ) ) {
			if ( new_e0.size() ) new_e0 = TrimCat( hbv.K(), new_e0, edges[ rpaths0.back() ] );
			else new_e0 = edges[ rpaths0.back() ];
			if ( new_e1.size() ) new_e1 = TrimCat( hbv.K(), new_e1, edges[ rpaths1.back() ] );
			else new_e1 = edges[ rpaths1.back() ];
			to_delete.push_back( from0 );
			to_delete.push_back( from1 );
			rpaths0.pop_back();
			rpaths1.pop_back();
        	    } else {
        		PRINT4( from0, from1, rpaths0.back(), rpaths1.back() );
        		FatalErr("BUG: lost synchromization with bubble edge list");
        	    }
        	}
                currv = hbv.From(currv)[0];
            }
            hbv.AddEdge(run.first, run.second, new_e0);
            hbv.AddEdge(run.first, run.second, new_e1);
            hbv.TestValid();
        }
    }
    hbv.DeleteEdges(to_delete);
    hbv.RemoveDeadEdgeObjects();
    hbv.RemoveEdgelessVertices();
}


int main(int argc, char *argv[])
{
    RunTime( );

    BeginCommandArguments;
    CommandDoc("Outputs edge adjacencies implied by the gap-free alignment of "
	    "a read pair to the HBV edges. Designed for a particular datatype "
	    "where the pairs imply potentially long-range relationships.");
    CommandArgument_String_Doc(QUERY_HEAD,
	    "specify the input paired reads QUERY_HEAD.{fastb,qualb,pairs}");
    CommandArgument_String_Doc(HBV, "binary file storing the HyperBaseVector");
    CommandArgument_String_OrDefault_Doc(OUT_HEAD,"" , "output OUT_HEAD.{hbv,dot}");
    CommandArgument_Int_OrDefault_Doc(L, 12, "flanking region used for hashing");
    CommandArgument_Bool_OrDefault(DEBUG, False);
    CommandArgument_Bool_OrDefault_Doc(REQUIRE_PERFECT, False, "require perfect matching reads");
    CommandArgument_Bool_OrDefault_Doc(PARTITION, False, "use partitioning algorithm instead");
    EndCommandArguments;

    std::cout << Date() << ": loading..." << endl;
    vecbasevector bases( QUERY_HEAD+".fastb");
    std::cout << Date() << ": loaded " << bases.size() << " reads." << endl;
    vecqualvector quals( QUERY_HEAD+".qualb");
    std::cout << Date() << ": loaded " << quals.size() << " quals." << endl;
    PairsManager pairs( QUERY_HEAD+".pairs");
    std::cout << Date() << ": PairsManager says " << pairs.nPairs() << " pairs and " << pairs.nReads() << " reads." << endl;
//    vecString names( QUERY_HEAD+".names" );

    HyperBasevector hbv;
    BinaryReader::readFile( HBV, &hbv );

    if ( OUT_HEAD != "" ) {
	Ofstream( dot, OUT_HEAD+".init.dot");
	hbv.PrintSummaryDOT0w(dot, true, true, true);
    }

//    VecAdjacencyVec adj = AlignPairs2( hbv, bases, quals, pairs, &names, REQUIRE_PERFECT, DEBUG );
    VecAdjacencyVec adj = AlignPairs2( hbv, bases, quals, pairs, 0, REQUIRE_PERFECT, DEBUG );

    if ( 1&&OUT_HEAD != "" ) {
	for ( size_t edge_no = 0; edge_no < adj.size(); ++edge_no ) {
	    for ( auto const& a : adj[edge_no] ) {
		cout << "EDGE " << edge_no << ": ";
		cout << a.edge2 << " [" << a.read1 << "," << a.read2 << "]";
		cout << endl;
	    }
	}
    }

    auto runs = find_bubble_runs(hbv, DEBUG);

    cout << "BUBBLE RUNS FOUND:" << endl;
    for ( auto const& run : runs ) cout << run.first << "," << run.second << endl;

    if ( PARTITION ) {
	vec<int> parts;
	part_bubble_runs_on_evidence( hbv, adj, runs, parts, DEBUG );

	cout << "BUBBLE RUNS FOUND AFTER PARTITIONING:" << endl;
	for ( auto const& run : runs ) cout << run.first << "," << run.second << endl;

	edit_bubble_runs_based_on_parts( hbv, runs, parts, DEBUG );
    }
    else {
	split_bubble_runs_on_evidence( hbv, adj, runs, DEBUG );

	cout << "BUBBLE RUNS FOUND AFTER SPLITTING:" << endl;
	for ( auto const& run : runs ) cout << run.first << "," << run.second << endl;

	edit_bubble_runs( hbv, runs, adj, DEBUG );
    }


    if ( OUT_HEAD != "" ) {
	Ofstream( dot, OUT_HEAD+".dot");
	hbv.PrintSummaryDOT0w(dot, true, true, true);
	hbv.DumpFasta(OUT_HEAD+".fasta",false);
    }
}

#include "feudal/SmallVecDefs.h"
template class SmallVec< Adjacency, MempoolAllocator<Adjacency> >;

#include "feudal/OuterVecDefs.h"
template class OuterVec<AdjacencyVec>;

template <> struct Serializability<Adjacency>
{ typedef TriviallySerializable type; };
