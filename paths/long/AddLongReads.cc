//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#include <omp.h>
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "paths/long/AddLongReads.h"

#include "Basevector.h"
#include "ParallelVecUtilities.h"
#include "paths/CorrectLongReads1.h"
#include "paths/CorrectLongReadsTools.h"
#include "paths/Uniseq.h"
#include "paths/GetNexts.h"
#include "paths/UnibaseUtils.h"
#include "paths/Uniseq.h"
#include "paths/HyperBasevector.h"
#include "paths/KmerBaseBroker.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/Unipath.h"
#include "paths/UnibaseUtils.h"
#include "math/Array.h"
#include "util/TextTable.h"
#include "kmers/KMer.h"
#include "paths/long/KmerAlign.h"
#include "paths/long/LongReadsToPaths.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "PrintAlignment.h"
#include "feudal/TrackingAllocator.h"
#include <queue>

namespace { // local utilities

template <int K, typename KmerT> 
void MakeKmerLookup(const vecbasevector& unibases, vec< triple<KmerT,int,int> >& kmers_plus)
{    vec<int64_t> starts;
     starts.push_back(0);
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    const basevector& u = unibases[i];
          starts.push_back( starts.back( ) + Max( 0, u.isize( ) - K + 1 ) );    }
     kmers_plus.resize( starts.back( ) );
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    const basevector& u = unibases[i];
          for ( int j = 0; j <= u.isize( ) - K; j++ )
          {    int64_t r = starts[i] + j;
               kmers_plus[r].first.SetToSubOf(u, j, K);
               kmers_plus[r].second = i; 
               kmers_plus[r].third = j;    }    }
     Sort(kmers_plus); 
}

HyperBasevector ReadsToHyperBaseVector(int K, const vecbasevector& reads,
                                        bool useOldLRPMethod )
{
    // generate unibases
    Bool ALT_PATHER= False;
    if ( ALT_PATHER ) {
        unsigned const COVERAGE = 2u;
        int verbosity = 0;
        HyperBasevector hbv;
        LongReadsToPaths(reads, K, COVERAGE, verbosity, useOldLRPMethod, &hbv);
        return hbv;
    }
    vecbasevector reads_mutable(reads);
    vecKmerPath paths, paths_rc;
    vec<tagged_rpint> pathsdb;
    ReadsToPathsCoreY( reads_mutable, K, paths, paths_rc, pathsdb, "",
                        getConfiguredNumThreads() );
    cout << Date( ) << ": calling Unipath" << endl;
    vecKmerPath unipaths;
    vec<tagged_rpint> unipathsdb;
    Unipath( paths, paths_rc, pathsdb, unipaths, unipathsdb );
    digraph A;
    cout << Date( ) << ": building adjacency graph" << endl;
    BuildUnipathAdjacencyGraph( paths, paths_rc, pathsdb, unipaths, unipathsdb, A );
    HyperKmerPath h;
    BuildUnipathAdjacencyHyperKmerPath( K, A, unipaths, h );
    KmerBaseBroker kbb( K, paths, paths_rc, pathsdb, reads );
    return HyperBasevector( h, kbb );
}

// Wrapper to Phase1 in CorrectLongReads
void PatchUsingPacbioReads(const int K, const vecbasevector& longreads, const vecbasevector& unibases,
                 vec<Bool> *p_computed, vec<GapPatcher> *p_patcher, vec<uniseq> *p_all_uniseqs,    
                 vec<vec<pair<int, vec<pair<int,int>>>>>* p_all_aligns,
                 int verbosity = 1 )
{
    const Bool PATCHE_ONLY = ( p_patcher != NULL );

    const int NUM_THREADS = 32;
    const int MIN_SAFE_CN1 = 150;
    const double MIN_VOTES_TO_PROTECT = 20;
    const Bool USE_SHORTEST = False;
    const Bool CLUSTER_ALIGNS_NEW_CLEAN = False;
    const Bool BEST_ONLY = False;
    const int MIN_PATCH1 = 200;
    const int MIN_PATCH2 = 200;
    const Bool STANDARD_ALIGNS = False;
    const Bool FILTER0 = False;
    const Bool FILTER = True;
    const Bool FILTER2 = False;
    const Bool CLEAN_GRAPH = False;
    const Bool CYCLIC_BAIL = True;
    const Bool KILL_INFERIOR_CLUSTERS_NEW = False;
    const Bool SCREEN_NEXTS = False;
    const int MIN_SPREAD1 = 10;
    const int MAX_GRAPH_SIZE = 0;
    const Bool SKIP_HIGHLY_COMPLEX_GRAPHS = True ;
    const Bool VERBOSE1 = ( verbosity >= 2 ? True: False ); // Default is False
    const Bool SKIP_SILENT = True;
    const Bool DOT1 = False;
    const Bool QLT1 = False;
    const Bool ABBREVIATE_ALIGNMENTS = True;
    const Bool PRINT_RAW_ALIGNS = False;
    const String data_dir = ""; // not used unless QLT1 == True in in Phase1() 
    const String run_dir = "";  // not used unless VERBOSE1y is set to >= 3 in Phase1()
    const Bool NEW_FILTER = False;
    heuristics heur;
    {
        heur.L = 11;
        heur.min_overlap_to_see_other = 250;
        heur.max_offset_diff_for_other = 25;
        heur.min_rdist_for_other = 10;
        heur.delta_r_max = 350;
        heur.max_delta_ratio = 1.5;
        heur.delta_sub = 10;
        heur.min_spread = 10;
        heur.min_ratio_to_kill = 5.0;
        heur.max_paths = 10000;
        heur.max_iterations = 10000;
        heur.min_safe_cn1 = MIN_SAFE_CN1;
        heur.min_win_ratio = 2.0;
        heur.min_votes = 4.0;
        heur.min_votes_to_protect = MIN_VOTES_TO_PROTECT;
    }
    const int L = heur.L;

    if ( verbosity >= 1 ) {
        cout << "K= " << K << endl;
        cout << "longreads.size()= " << longreads.size() << endl;
        cout << "nuni= " << unibases.size() << endl;
    }
    vec<int> reads_to_use(longreads.size(), vec<int>::IDENTITY); 
    vec<int> reads_to_trace; // empty: do not trace

    cout << Date( ) << ": making nexts" << endl;
    vec< vec<int> > nexts;
    GetNexts( K, unibases, nexts );
    vec<int> to_rc;
    UnibaseInvolution( unibases, to_rc );

    int nuni = unibases.size( ); 
    vec< vec<int> > fillseqs;
    vec< vec<int> > nexts_count(nuni);
    {
        cout << Date() << ": creating ffpairs" << endl;
        vec< pair<int,int> > ffpairs;
        for ( int i = 0; i < fillseqs.isize( ); i++ )
        {    const vec<int>& x = fillseqs[i];
            for ( int j = 0; j < x.isize( ) - 1; j++ )
                ffpairs.push( x[j], x[j+1] );    }
        cout << Date( ) << ": sorting" << endl;
        ParallelSort(ffpairs);
        cout << Date() << ": cataloging" << endl;
        for ( int u = 0; u < nuni; u++ )
            nexts_count[u].resize( nexts[u].size( ), 0 );
        for ( int i = 0; i < ffpairs.isize( ); i++ )
        {    int j = ffpairs.NextDiff(i);
            int u = ffpairs[i].first, v = ffpairs[i].second;
            for ( int l = 0; l < nexts[u].isize( ); l++ )
                if ( nexts[u][l] == v ) nexts_count[u][l] = j - i;
            i = j - 1;    }
     }

    cout << Date() << ": hashing unibases" << endl;
    vec< vec< pair<int,int> > > Ulocs( IPow( 4, L ) );
    for ( size_t i = 0; i < unibases.size( ); i++ ) 
    {    for ( int j = 0; j <= unibases[i].isize( ) - L; j++ ) 
          {    int n = KmerId( unibases[i], L, j );
               Ulocs[n].push( i, j );    }    } 

    int nreads = longreads.size( );
    vec<uniseq> UNISEQ(nreads);
    vec< vec<int> > UNISEQ_ID(nreads);
    vec< vec< vec<int> > > BESTS(nreads);
    vec<Bool> COMPUTED( nreads, False );
    vecbasevector all;
    vec<GapPatcher> patchers;
    cout << Date() << ": START PHASE 1, PASS 1" << endl;
    vec< digraphVE<int,int> > Hall(nreads);
    vec< vec< pair< int, vec< pair<int,int> > > > > ALIGNS_ALL(nreads);
    Phase1( NUM_THREADS, unibases, to_rc, nexts, nexts_count, K, L, Ulocs, 
          longreads, reads_to_use, heur, USE_SHORTEST, CLUSTER_ALIGNS_NEW_CLEAN,
          BEST_ONLY, PATCHE_ONLY, MIN_PATCH1, MIN_PATCH2, STANDARD_ALIGNS, FILTER0, 
          FILTER, FILTER2, CLEAN_GRAPH, CYCLIC_BAIL, KILL_INFERIOR_CLUSTERS_NEW, 
          SCREEN_NEXTS, MIN_SPREAD1, MAX_GRAPH_SIZE, SKIP_HIGHLY_COMPLEX_GRAPHS, reads_to_trace, 
          VERBOSE1, SKIP_SILENT, DOT1, QLT1, data_dir, run_dir, 
          ABBREVIATE_ALIGNMENTS, PRINT_RAW_ALIGNS, UNISEQ, UNISEQ_ID, BESTS, 
          COMPUTED, all, patchers, Hall, ALIGNS_ALL );

    *p_computed = COMPUTED;
    *p_all_uniseqs = UNISEQ;
    *p_all_aligns = ALIGNS_ALL;
    if ( p_patcher != NULL ) *p_patcher = patchers;
}

// Find all paths (upto max_num_paths) in the graph from left_edge to
// right_dege closest to the expected gap size. Do not try path larger than
// upper_bound. Quit search if more than 10000 paths has been tested.
vec<vec<int>> GetPathsForGap( const HyperBasevector& shb, int left_edge, int right_edge, 
        int gap_expect, int upper_bound, int max_num_paths = 5000) 
{
    const int max_to_generate = 10000;
    priority_queue<pair<int,vec<int>>> top_paths;
    const int KS = shb.K();
    vec<int> to_right;
    shb.ToRight(to_right);
    vec<int> to_left;
    shb.ToLeft(to_left);
    // The search will stop if the path pass the right_edge
    int v = to_right[ left_edge ];
    int w = to_left[ right_edge ];
    vec<vec<int>> partials;
    partials.push( vec<int>() );
    int ntested=0;
    while( ! partials.empty() ) {    
        //if ((int)top_paths.size() >= max_num_paths) break;
        if (ntested > max_to_generate) break;
        vec<int> p = partials.back();
        partials.pop_back( );
        int gap_size = -KS+1;
        for ( int i = 0; i < p.isize( ); i++ ) {
            gap_size += shb.EdgeLength(p[i]) - KS + 1;
        }
        int gap_dev = abs(gap_size - gap_expect);
        int vn = ( p.empty( ) ? v : to_right[ p.back( ) ] );
        if ( vn == w ) {
            ntested++;
            if ((int)top_paths.size() < max_num_paths || top_paths.top().first < gap_dev){
                if ((int)top_paths.size() >= max_num_paths) 
                    top_paths.pop();
                top_paths.push(make_pair(gap_dev, p)); // destination node reached   
            }
        }
        for ( int j = 0; j < shb.From(vn).isize( ); j++ ) {    
            int e = shb.EdgeObjectIndexByIndexFrom( vn, j );
            if (gap_size + shb.EdgeLength(e) - KS + 1 > upper_bound ) 
                continue;
            vec<int> p2 = p;
            p2.push_back(e);
            partials.push_back(p2);    
        } 
    }
    vec<vec<int>> paths(top_paths.size());
    for (auto it = paths.rbegin(); it != paths.rend(); ++it) {
       *it = top_paths.top().second; 
       top_paths.pop();
    }
    return paths;
}
// Exaustively find all edges that can be visited starting from left_edge without
// passing right_ledge. 
vec<int> GetEdgesBetween( const SupportedHyperBasevector& shb, int left_edge, int right_edge)
{
    vec<int> to_right;
    shb.ToRight(to_right);
    // The search will stop if the path pass the right_edge
    vec<Bool> visited(shb.EdgeObjectCount(), False);
    vec<int> stack = { left_edge };
    while( ! stack.empty() ) {    
        int edge = stack.back();
        stack.pop_back();
        visited[edge] = True;
        if (edge == right_edge) continue; // do not pass right_edge
        int vn = to_right[edge];
        for (int j = 0; j < shb.From(vn).isize( ); j++) {    
            int e = shb.EdgeObjectIndexByIndexFrom(vn, j);
            if (!visited[e])
                stack.push_back(e);
        } 
    }
    vec<int> vec_edge;
    for (int i = 0; i < visited.isize(); ++i) 
        if (visited[i] && i != left_edge && i != right_edge) 
            vec_edge.push_back(i);
    return vec_edge;
}

vec<int> GetEdgesDownStream( const SupportedHyperBasevector& shb, int left_edge, int right_edge)
{
    vec<int> to_right;
    shb.ToRight(to_right);
    // The search will stop if the path pass the right_edge
    vec<Bool> visited(shb.EdgeObjectCount(), False);
    vec<int> stack = { left_edge };
    while( ! stack.empty() ) {    
        int edge = stack.back();
        stack.pop_back();
        visited[edge] = True;
        if (edge == right_edge) {
            visited[right_edge] = True;
            continue; 
        }
        int vn = to_right[edge];
        for (int j = 0; j < shb.From(vn).isize( ); j++) {    
            int e = shb.EdgeObjectIndexByIndexFrom(vn, j);
            if (!visited[e])
                stack.push_back(e);
        } 
    }
    vec<int> vec_edge;
    for (int i = 0; i < visited.isize(); ++i) 
        if (visited[i]) vec_edge.push_back(i);
    Sort(vec_edge);
    return vec_edge;
}

vec<int> GetEdgesUpStream( const SupportedHyperBasevector& shb, int right_edge, int left_edge)
{
    vec<int> to_left;
    shb.ToLeft(to_left);
    vec<Bool> visited(shb.EdgeObjectCount(), False);
    vec<int> stack = { right_edge };
    while( ! stack.empty() ) {    
        int edge = stack.back();
        stack.pop_back();
        visited[edge] = True;
        if (edge == left_edge) {
            visited[left_edge] = True;
            continue; 
        }
        int vn = to_left[edge];
        for (int j = 0; j < shb.To(vn).isize( ); j++) {    
            int e = shb.EdgeObjectIndexByIndexTo(vn, j);
            if (!visited[e])
                stack.push_back(e);
        } 
    }
    vec<int> vec_edge;
    for (int i = 0; i < visited.isize(); ++i) 
        if (visited[i]) vec_edge.push_back(i);
    Sort(vec_edge);
    return vec_edge;
}

// Find the cells that we want to unroll from the graph.
// We want to unroll the edges between two edges e1 and e2 if:
// 1. Both are long edges (defined as >= 1K bases or >=150 for end edges)
// 2. There are multiple paths go through these two edge. 
//      path1: e1-> ... -> e2
//      path2: e1-> ... -> e2
// 3. All paths from e1 always gos to e2 and vice versa
vec<pair<int,int>> FindUnrollingPlaces(const SupportedHyperBasevector& shb)
{
    int verbosity = 1;
    const int LongEdgeLenMin = 300;
    const int EndEdgeLenMin = 150;
    vec<int> to_right;  
    shb.ToRight(to_right);
    vec<int> to_left;
    shb.ToLeft(to_left);
    // find end edges
    vec<int> endEdges;
    for (int i = 0; i < shb.N(); ++i) 
        if (shb.To(i).size() == 1 && shb.From(i).size() == 0)
            endEdges.push_back(shb.EdgeObjectIndexByIndexTo(i,0));
        else if (shb.From(i).size() == 1 && shb.To(i).size() == 0)
            endEdges.push_back(shb.EdgeObjectIndexByIndexFrom(i,0));
    UniqueSort(endEdges);
    endEdges.Println(cout);
    vec<int> longEdgeVec;
    for ( int i = 0; i < shb.EdgeObjectCount( ); i++ ) {
        if (shb.EdgeLength(i) >= LongEdgeLenMin ||
             ((int)shb.EdgeLength(i) >= EndEdgeLenMin 
               && BinMember(endEdges,i)))
            longEdgeVec.push_back(i);
    }
    longEdgeVec.Println(cout);
    UniqueSort(longEdgeVec);
    // Find unrolling edge pairs from longEdgeVec
    vec<pair<int,int>> edgePairs;
    vec<vec<int>> inner_edges;
    for( int edge1 : longEdgeVec) {
    for( int edge2 : longEdgeVec) {
        if (edge1 == edge2) continue;
        //if (to_right[edge1] == to_left[edge2]) continue; // edge1-> x -> edge2
        // If downstream from edge1 and upstream from edge2 are the same
        vec<int> down_stream = GetEdgesDownStream(shb, edge1, edge2);
        vec<int> up_stream = GetEdgesUpStream(shb, edge2, edge1);
        if (down_stream == up_stream) {
            edgePairs.push(edge1,edge2);
            inner_edges.push_back(down_stream);
        }
    }}
    for (auto x: edgePairs)
        cout << "candidate path: " << x.first << "," << x.second << endl;
    vec<Bool> todel(edgePairs.size(),False);
    for (size_t i = 0; i < inner_edges.size(); ++i) {
    for (size_t j = 0; j < inner_edges.size(); ++j) {
        if (j == i) continue;
        if ( Member(inner_edges[i],edgePairs[j].first) &&
             Member(inner_edges[i],edgePairs[j].second) ) {
            todel[i] = True;
            if (verbosity>=1)
                cout << "path between " << edgePairs[i].first << "," << edgePairs[i].second
                    << " contains " << edgePairs[j].first << "," << edgePairs[i].second << endl;
            break;
        }
    }}
    EraseIf(edgePairs,todel);
    EraseIf(inner_edges,todel);
    if (verbosity >=1 ) {
        cout << "finding " << edgePairs.size() << " pairs: " << endl;
        for(auto x: edgePairs) cout << x.first << " " << x.second << endl;
    }
    return edgePairs;
}

// Unrolling by forming pacbio consensus read
basevector PacbioConsensus(int K, basevector left_seq, basevector right_seq, 
        const vecbasevector& longreads, bool Trim_Repeat_Ends = false)
{
    const int verbosity = 2;
    // trim off any highly repetitive ends
    if ( Trim_Repeat_Ends ) {
        vecbasevector two_reads;
        two_reads.push_back(left_seq);
        two_reads.push_back(right_seq);
        const int KS = 8;
        vec<triple<basevector,int,int>> kmers_plus;
        MakeKmerLookup<KS,basevector>(two_reads, kmers_plus);
        vec<vec<bool>> repeat_pos(2);
        repeat_pos[0].resize(left_seq.size() - KS + 1, false);
        repeat_pos[1].resize(right_seq.size() - KS + 1, false);
        for (size_t i = 0; i < kmers_plus.size(); ++i) {
            size_t j = i+1;
            while (j < kmers_plus.size() && kmers_plus[j].first == kmers_plus[i].first)
                j++;
            if (j - i > 1) 
                repeat_pos[kmers_plus[i].second][kmers_plus[i].third] = true;
        }
        int left_repeat_trim = 0, right_repeat_trim = 0;
        for (auto it = repeat_pos[0].rbegin(); it != repeat_pos[0].rend(); ++it) 
            if (*it) left_repeat_trim++;
        for (auto it = repeat_pos[1].begin(); it != repeat_pos[1].end(); ++it) 
            if (*it) right_repeat_trim++;
        cout << "left_repeat_trim= " << left_repeat_trim << endl;
        cout << "right_repeat_trim= " << right_repeat_trim << endl;
        left_seq.SetToSubOf(left_seq, 0, left_seq.size() - left_repeat_trim);
        right_seq.SetToSubOf(right_seq, right_repeat_trim, right_seq.size() - right_repeat_trim);
    }
    //Patch using pacbio reads
    basevector result;
    {
        vecbasevector edge2base_local;
        edge2base_local.push_back(left_seq);
        edge2base_local.push_back(right_seq);
        vecbasevector unibases;
        //ReadsToUnibases( K, edge2base_local, &unibases );
        unibases.push_back(left_seq);
        unibases.push_back(right_seq);
        unibases.push_back(ReverseComplement(left_seq));
        unibases.push_back(ReverseComplement(right_seq));
        vec<Bool> computed;
        vec<GapPatcher> patchers;
        vec<uniseq> all_uniseqs;
        vec<vec<pair<int, vec<pair<int,int>>>>> all_aligns;
        // find all patches and keep only patches bridging the two edges
        PatchUsingPacbioReads(K, longreads, unibases, &computed, &patchers, &all_uniseqs, &all_aligns,1);
        if (verbosity >= 1) {
            cout << "Had " << computed.CountValue(True) << " computed" << endl;
            cout << "Found " << patchers.size() << " patches" << endl;
            //map<pair<int,int>,int> patcher_stats;
            //for (size_t i = 0; i < patchers.size(); ++i) {
            //    patcher_stats[make_pair(patchers[i].t1, patchers[i].t2)]++;
            //}
            //for(auto x: patcher_stats)
            //    cout << x.first.first << " " << x.first.second << " " << x.second << endl;
        }
        vec<Bool> dels(patchers.size(), True);
        for (size_t i = 0; i < patchers.size(); ++i) {
            if(patchers[i].t1 == 0 && patchers[i].t2 == 1) dels[i] = False;
        }
        EraseIf(patchers, dels);
        if (verbosity >= 1) {
            for (size_t i = 0; i < patchers.size(); ++i) {
                // the ids saved in GapPatcher are indexed with fw + rc interlaced
                int pacbio_rid = patchers[i].rid / 2;
                int rc = patchers[i].rid % 2;
                if (verbosity >= 2)
                    cout << pacbio_rid << (rc? " ": " ") << " rsize=" <<
                        longreads[pacbio_rid].size() << " tpos=" << patchers[i].tpos1 << " "
                        << patchers[i].tpos2 << " gap=" << patchers[i].expected_gap << endl;
            }
        }
        // Build consensus using the patches
        Bool CORRECT_PATCHES = False;
        Bool CORRECT_PATCHES_NEW = False;
        Bool CORRECT_PATCHES_VERBOSE= False;
        Bool ANNOUNCE_PATCH_CORRECTION= False;
        vec<Bool> cn_plus(unibases.size(), False);  // force copy number to be one
        const int L = 11;
        const int LG = 12;
        cout << Date( ) << ": hashing unibases" << endl;
        vec< vec< pair<int,int> > > Ulocs( IPow( 4, L ) );
        for ( size_t i = 0; i < unibases.size( ); i++ ) 
        {    for ( int j = 0; j <= unibases[i].isize( ) - L; j++ ) 
            {    int n = KmerId( unibases[i], L, j );
                 Ulocs[n].push( i, j );    }    } 
        // Below are needed only if validation is required
        vecbasevector genome2;
        vec< vec< pair<int,int> > > Glocs;
        Bool PATCH_VERBOSITY = 0;
        Bool VALIDATE_PATCHES = False;
        Bool PATCH_CORRECT_VERBOSITY= 0;
        Bool MIN_TO_PATCH= 5;
        String data_dir = "";
        String run_dir = "";
        vecbasevector fbases;
        vecqualvector fquals;
        PairsManager fpairs;

        const int F = 20;
        cout << Date( ) << ": building index" << endl;
        vec< kmer<F> > fheads( fbases.size( ) );
        for ( size_t id = 0; id < fbases.size( ); id++ )
            if ( fbases[id].isize( ) >= F ) fheads[id].SetToSubOf( fbases[id], 0 );
        vec<int64_t> fids( fbases.size( ), vec<int64_t>::IDENTITY );
        ParallelSortSync( fheads, fids );
        Bool NEW_PLUS_RULE= False;

        cout << Date( ) << ": making nexts" << endl;
        vec< vec<int> > nexts;
        GetNexts( K, unibases, nexts );
        vec<basevector> bpatches;
        cout << Date() << ": Building Patches <96>" << endl;   
        BuildPatches<96>( unibases, nexts, cn_plus, L, Ulocs, patchers, 
                CORRECT_PATCHES, CORRECT_PATCHES_NEW,
                CORRECT_PATCHES_VERBOSE, ANNOUNCE_PATCH_CORRECTION, genome2, 
                PATCH_VERBOSITY, VALIDATE_PATCHES, PATCH_CORRECT_VERBOSITY, LG, 
                Glocs, data_dir, run_dir, bpatches, MIN_TO_PATCH,
                fbases, fquals, fpairs, fheads, fids, NEW_PLUS_RULE );    
        cout << Date() << ": Patches built <96>" << endl;
        if (verbosity >= 1) {
            cout << "left_seq= " << left_seq << endl;
            cout << "right_seq= " << right_seq << endl;
            for (auto x: bpatches)
                cout << "patch= " << x << endl;
        }
        if ( bpatches.size() == 1 )
            result = bpatches[0];
    }
    return result;
}

// Join continuous edges
// v1 -> v2 -> v3 -> v4 ====> v1 -> v4
int JoinConsecutiveEdges(SupportedHyperBasevector& shb)
{
    int verbosity = 1;
    int njoin = 0;
    while (true) {
        bool changed = false;
        for(int v = 0; v < shb.N(); v++) {
            if (shb.From(v).solo() && shb.To(v).solo() ) {
                int e1 = shb.EdgeObjectIndexByIndexTo(v, 0 );
                int e2 = shb.EdgeObjectIndexByIndexFrom(v, 0 );
                basevector base1 = shb.EdgeObject(e1);
                basevector base2 = shb.EdgeObject(e2);
                base1.resize(base1.size() - shb.K() + 1);
                shb.JoinEdges(v, Cat(base1, base2));
                changed = true;
                if (verbosity >= 1)
                    cout<< "join edges " << e1 << " and " << e2 << " at node " << v << endl;
                njoin++;
                break;
            }
        }
        if ( !changed) break;
    }
    return njoin;
}

template<typename GraphT>
void DumpGraph(const String filename, const GraphT& shb) {
    // dump the simplified graph
    vec<String> edge_names2(shb.EdgeObjectCount(),"");
    vec<double> lengths2( shb.EdgeObjectCount( ) );
    for (size_t i = 0; i < edge_names2.size(); ++i) {
        edge_names2[i] = ToString(i) + "len= " + ToString(shb.EdgeObject(i).size());
        lengths2[i] = shb.EdgeLengthKmers(i);
    }
    ofstream dout(filename);
    shb.PrettyDOT( dout, lengths2, HyperBasevector::edge_label_info(
                HyperBasevector::edge_label_info::DIRECT, &edge_names2 ) );
}


// Find the sequence and matching score from a graph that matches best the given read.
pair<int, basevector> BestSeqFromGraph(const basevector& bb, const HyperBasevector& hb, 
         int left_edge, int left_edge_trim, int right_edge, int right_edge_trim,
         int max_errors)
{
    int verbosity = 1;
    const int max_num_paths = 5000;
    int upper_bound = bb.size() + max_errors;
    int lower_bound = bb.size() - max_errors;
    const int KS = hb.K();
    int gap_size = bb.isize() - hb.EdgeLength(left_edge) - hb.EdgeLength(right_edge)
        + left_edge_trim + right_edge_trim;
    cout << "gap_size= " << gap_size << endl;
    vec<vec<int>> all_paths = GetPathsForGap( hb, left_edge, right_edge, 
            gap_size, gap_size + max_errors, max_num_paths);
    if (verbosity >= 1 )
        cout << "number of paths prepared " << all_paths.size() << endl;;

    vec<basevector> reads_from_paths(all_paths.size());
    #pragma omp parallel for
    for(size_t i = 0; i < all_paths.size(); i++) {
        vec<int> path = all_paths[i];
        path.insert(path.begin(), left_edge);
        path.insert(path.end(), right_edge);
        basevector x = (path.empty()? basevector() : hb.EdgePathToBases(path));
        basevector x2(x.begin() + left_edge_trim, x.end() - right_edge_trim);
        //ForceAssertLe(lower_bound, x2.isize());
        //ForceAssertLe(x2.isize(), upper_bound);
        reads_from_paths[i] = x2;
    }
    if (verbosity >= 1 )
        cout << "number of reads to be tested " << reads_from_paths.size() << endl;;
    vec<int> scores(reads_from_paths.size(), -1);
    bool UseRecursiveAligner = false;
    if (UseRecursiveAligner) {
        int I_Score = 1, D_Score = 1, S_Score = 1, U_Score = 0;
        int MaxSeqLen = 0;
        for (const auto &x: reads_from_paths) MaxSeqLen = Max(MaxSeqLen, (int)x.size());
        int M = MaxSeqLen, N = bb.size();
        RecArray<int> s( M+1, N+1 ); 
        cout << Date() << ": Aligning read to graph." << endl;
        for (int i = 0; i <= M; ++i) s[i][0] = D_Score * i;
        for (int j = 0; j <= N; ++j) s[0][j] = I_Score * j;
        int best_score_track = 10000000;
        int break_pos = M;
        for (size_t rid = 0; rid < reads_from_paths.size(); ++rid) {
            const basevector& aa = reads_from_paths[rid];
            if (verbosity >= 1)
                if (rid % (reads_from_paths.isize()/100) == 0 ) cout << '.' << flush;
            if (verbosity >= 2 )
                cout << "rid= " << rid << " aa.size()= " << aa.size() << " bb.size()= " << bb.size() << endl;
            int last_diff_pos = 0;
            if (rid>0) {
                size_t x1 = 0, x2 = 0;
                while ( x1 != aa.size() && x2 != reads_from_paths[rid-1].size() && 
                        aa[x1] == reads_from_paths[rid-1][x2] )
                    ++last_diff_pos, ++x1, ++x2;
                last_diff_pos = min(last_diff_pos, break_pos);
            }
            bool bad_path = false;
            for (int i = last_diff_pos+1; i <= (int)aa.size(); i++) {
                int min_score = s[i][0];
                for (int j = 1; j <= N; j++) {
                    s[i][j] = std::min({aa[i-1] == bb[j-1] ? s[i-1][j-1] + U_Score: s[i-1][j-1] + S_Score,
                            s[i-1][j] + D_Score,  s[i][j-1] + I_Score});
                    min_score = min(min_score, s[i][j]);
                }
                break_pos = i;
                if (min_score > best_score_track) {
                    s[aa.size()][N] = 10000000;
                    bad_path = true;
                    break;
                }
            }
            if (bad_path)
                scores[rid] = 10000000;
            else {
                scores[rid] = s[aa.size()][N];
                best_score_track = scores[rid];
            }
        }
        if (verbosity >= 1) cout << endl;
    } else {
        #pragma omp parallel for
        for (size_t rid = 0; rid < reads_from_paths.size(); ++rid) {
            align a;
            SmithWatFreeSym(reads_from_paths[rid], bb, a, true, true, 1, 1, 0);
            scores[rid] = a.Errors(reads_from_paths[rid], bb);
        }
    }
    cout << Date() << ": Looking for best match." << endl;
    int best_score = scores[0], best_score_index = 0;
    for (size_t i = 0; i < scores.size(); ++i) {
        if (scores[i] < best_score) {
            best_score = scores[i];
            best_score_index = i;
        }
    }
    if (verbosity >= 1 )
        cout << "best_score= " << best_score << " best_score_index= " << best_score_index << endl;
    return make_pair(best_score, reads_from_paths[best_score_index]);
}

// Create the local graph from the edges and locate the position of the starting kmer of
// the left_seq and ending kmer of the right_seq.
void CreateLocalGraphAndLocate(const int KS, const SupportedHyperBasevector& shb, 
        const basevector& left_seq, const basevector& right_seq, const vec<int>& edges,
        HyperBasevector* jungle_hbv, int* left_edge2, int* left_edge2_trim,
        int* right_edge2, int* right_edge2_trim, bool useOldLRPMethod )
{
    basevector head_kmer(left_seq.begin(), left_seq.begin()+KS);
    basevector tail_kmer(right_seq.end() - KS, right_seq.end());
    ForceAssertNe(head_kmer, tail_kmer); // don't handle this case yet
    vecbasevector jungle;
    for (auto x: edges) jungle.push_back(shb.EdgeObject(x));
    jungle.push_back(left_seq);
    jungle.push_back(right_seq);
    *jungle_hbv = ReadsToHyperBaseVector(KS, jungle, useOldLRPMethod);
    *left_edge2 = *left_edge2_trim = *right_edge2 = *right_edge2_trim = -1;
    for (int i = 0; i < jungle_hbv->EdgeObjectCount(); ++i) {
        for (int j = 0; j < jungle_hbv->EdgeLength(i) - KS + 1; ++j) {
            basevector e1(jungle_hbv->EdgeObject(i), j, KS);
            if (e1 == head_kmer) {
                *left_edge2 = i;
                *left_edge2_trim = j;
                break;
            }
        }
        if (*left_edge2 != -1) break;
    }
    for (int i = 0; i < jungle_hbv->EdgeObjectCount(); ++i) {
        for (int j = 0; j < jungle_hbv->EdgeLength(i) - KS + 1; ++j) {
            basevector e1(jungle_hbv->EdgeObject(i), jungle_hbv->EdgeLength(i) - KS -j, KS);
            if (e1 == tail_kmer) {
                *right_edge2 = i;
                *right_edge2_trim = j;
                break;
            }
        }
        if (*right_edge2 != -1) break;
    }
}

vec<int> FindLongSingletonEdges(const SupportedHyperBasevector& shb, const int LongEdgeLenMin = 1000)
{
    int verbosity = 1;
    vec<int> to_right;  
    shb.ToRight(to_right);
    vec<int> to_left;
    shb.ToLeft(to_left);
    vec<int> longEdgeVec;
    for ( int i = 0; i < shb.EdgeObjectCount( ); i++ )
        if (shb.EdgeLength(i) >= LongEdgeLenMin &&
                shb.From(to_left[i]).solo() && shb.To(to_left[i]).empty() &&
                shb.From(to_right[i]).empty()&&shb.To(to_right[i]).solo() )
            longEdgeVec.push_back(i);
    if (verbosity >=1 )
        longEdgeVec.Println(cout);
    return longEdgeVec;
}

} // end anonymous namespace

int PatchWithPacbioReads(const vecbasevector& longreads, SupportedHyperBasevector * p_shb) 
{
    int njoin = 0;
    const int verbosity = 2;
    SupportedHyperBasevector& shb = *p_shb;
    const int K = shb.K();
    vec<int> candidate_edges = FindLongSingletonEdges(shb);
    typedef std::basic_ostringstream<char,std::char_traits<char>,typename DefaultAllocator<char>::type> StdOSS;
    vec<StdOSS> status_logs(candidate_edges.size());
    const int edge_trim_len = 500; // only keep max 500 base edge for anchoring pacbio reads
    for(size_t i = 0; i < candidate_edges.size(); i++) {
        if (verbosity >= 1)
            cout << "================= Try to extend edge ================ " << i << endl;
        int left_edge = candidate_edges[i];
        int right_edge = -1, ltrim = -1, rtrim = -1;
        basevector consensus, ll, rr, left_seq, right_seq;
        ll = shb.EdgeObject(left_edge); 
        //bool rc = false;
        for(int pass = 0; pass < 2; pass++) {
            if (pass == 1) 
                ll.ReverseComplement();
            for (auto x: candidate_edges) {
                if (x == left_edge) continue;
                rr = shb.EdgeObject(x);
                ltrim = min(edge_trim_len, ll.isize()); // actual trimmed size
                rtrim = min(edge_trim_len, rr.isize());
                left_seq = basevector(ll.end() - ltrim, ll.end());
                right_seq = basevector(rr.begin(), rr.begin() + rtrim);
                consensus = PacbioConsensus(K, left_seq, right_seq, longreads);
                if (!consensus.empty()) {
                    //rc = (pass == 1);
                    right_edge = x;
                    break;
                }
            }
            if (! consensus.empty()) break;
        }
        if (consensus.empty()) {
            cout << "Unable to generate consensus" << endl;
            continue;
        }
        int matching_left = 0;
        for(int x = 0; x < min(ltrim, consensus.isize()); x++)
            if (left_seq[x] == consensus[x]) 
                matching_left++;
        int matching_right = 0;
        for(int x = 0; x < min(rtrim, consensus.isize()); x++)
            if (right_seq[rtrim - 1 - x] == consensus[consensus.size() - 1- x]) 
                matching_right++;
        BaseVec new_edge(consensus.begin() + matching_left - K + 1, 
                consensus.end() - matching_right + K - 1); // K-1 overlap with neighbor edges
        if (verbosity >=1) {
            cout << "matching_left= " << matching_left << " ltrim= " << ltrim
                << " matching_right= " << matching_right << " rtrim= " << rtrim 
                << " edge_len= " << new_edge.size() << endl;
        }
        vec<int> to_right; shb.ToRight(to_right);
        vec<int> to_left; shb.ToLeft(to_left);
        shb.EdgeObjectMutable(left_edge) = BaseVec(ll.begin(), ll.end() - ltrim + matching_left);
        shb.EdgeObjectMutable(right_edge) =
            basevector(rr.begin() + rtrim - matching_right, rr.end());
        shb.AddEdge(to_right[left_edge], to_left[right_edge], new_edge);
        cout << "New edge added, size= " << new_edge.size() << endl;
        njoin += JoinConsecutiveEdges(shb);
    }
    return njoin;
}

void UnrollWithPacbioReads(const vecbasevector& longreads, SupportedHyperBasevector * p_shb,
        bool useOldLRPMethod, bool CorrectPacbioConsensus, int verbosity)
{
    SupportedHyperBasevector& shb = *p_shb;
    if (verbosity >= 2) {
        BinaryWriter::writeFile("input.shbv", shb );
        BinaryWriter::writeFile("pacbio.fastb",longreads);
        DumpGraph("before3.dot", shb);
    }
    HyperBasevector hb(shb);
    const int K = hb.K();
    vec<int> to_right;
    shb.ToRight(to_right);
    vec<int> to_left;
    shb.ToLeft(to_left);

    vec<pair<int,int>> edge_pairs = FindUnrollingPlaces(shb);
    typedef std::basic_ostringstream<char,std::char_traits<char>,typename DefaultAllocator<char>::type> StdOSS;
    vec<StdOSS> status_logs(edge_pairs.size());
    const int edge_trim_len = 500; // only keep max 500 base edge for anchoring pacbio reads
    for(size_t ipair = 0; ipair < edge_pairs.size(); ipair++) {
        int left_edge = edge_pairs[ipair].first;
        int right_edge = edge_pairs[ipair].second;
        ostream & status_log = status_logs[ipair];
        cout << " ====================== Try to unroll between edge " << left_edge << 
            " and " << right_edge << " ========================   " << endl;
        basevector ll = shb.EdgeObject(left_edge), rr = shb.EdgeObject(right_edge);
        int ltrim = min(edge_trim_len, ll.isize()); // actual trimmed size
        int rtrim = min(edge_trim_len, rr.isize());
        basevector left_seq(ll.end() - ltrim, ll.end());
        basevector right_seq(rr.begin(), rr.begin() + rtrim);
        vec<int> edges = GetEdgesBetween(shb, left_edge, right_edge);
        bool has_long_edge = false;
        for (auto x: edges) 
            if (shb.EdgeLength(x) > 2000) {
                has_long_edge = true;
                break;
            }
        if (has_long_edge) {
            status_log << "Don't unroll  long edges" << endl;
            continue;
        }
        basevector consensus = PacbioConsensus(K, left_seq, right_seq, longreads);
        if (consensus.empty()) {
            status_log << "Unable to generate consussus";
            continue;
        }

        // If unrolling failed, don't do any editing. Otherwise, insert the
        // unrolled sequence between left and right edges. Note that the
        // unrolling may change the left and right edges. If that happens, the
        // two edges will be trimmed so that the changes will only be in the inserted edge.

        if (CorrectPacbioConsensus) {
            //// now test the consensus
            // ------------->* edges between*------------ right_seq
            // ------------------------------------------ consensus
            double max_error_rate = 0.05;
            int max_errors =  max(K, consensus.isize() - ltrim - rtrim) * max_error_rate;
            max_errors = max(max_errors, 5);
            vec<pair<int,basevector>> best_reads(2);
            // =============== First round ================
            if (verbosity >= 1) cout << "Validating consensus using original graph " << endl;
            int left_edge_trim = ll.isize() - left_seq.isize();
            int right_edge_trim = rr.isize() - right_seq.isize();
            best_reads[0] = BestSeqFromGraph(consensus, shb, left_edge, left_edge_trim, 
                    right_edge, right_edge_trim, max_errors); 
            // =============== Second round ================
            const int KS = 24;
            int left_edge2 = -1, left_edge2_trim = -1;
            int right_edge2 = -1, right_edge2_trim = -1;
            HyperBasevector jungle_hbv;
            CreateLocalGraphAndLocate(KS, shb, left_seq, right_seq, edges,
                    &jungle_hbv, &left_edge2, &left_edge2_trim, 
                    &right_edge2, &right_edge2_trim, useOldLRPMethod);

            if (verbosity >= 2)
                DumpGraph(ToString(left_edge2) + "-" + ToString(right_edge2) + "jungle.dot", jungle_hbv);
            if (verbosity >= 1) {
                cout << "Validating consensus using K=24 graph" << endl;
                cout << "Start from edge " << left_edge2 << " skip_head= " << left_edge2_trim 
                    << " to edge " << right_edge2 << " skip_tail= " << right_edge2_trim << endl;
            }
            best_reads[1] = BestSeqFromGraph(consensus, jungle_hbv, left_edge2, left_edge2_trim, 
                    right_edge2, right_edge2_trim, max_errors); 
            if (verbosity >= 1) {
                cout << "mismatches= " << best_reads[0].first << " " << best_reads[1].first << endl;
            }
            int best_pass = (best_reads[0].first < best_reads[1].first ? 0 : 1);
            basevector best_read_from_graph = best_reads[best_pass].second;
            align a;
            SmithWatFreeSym(best_read_from_graph, consensus, a, true, true, 1, 1, 0 );
            int nerror = a.Errors(best_read_from_graph, consensus);
            float error_rate = nerror * 1.0 / max(K, consensus.isize() - ltrim - rtrim);
            if (verbosity >=1 )
                cout << "error_rate= " << error_rate << " nerror= " << nerror << endl;
            PrintVisualAlignment(true, cout, best_read_from_graph, consensus, a);
            if (nerror > max_errors) {
                if (verbosity >= 1)
                    cout << "Reject consensus read nerror= " << nerror << " error_rate= " << error_rate << endl;;
                status_log << "Reject consensus read nerror= " << nerror << " error_rate= " << error_rate << endl;;
                continue;
            } else {
                if (verbosity >= 1){
                    cout << "Accept consensus read after making " << nerror << " corrections" << endl;
                }
                if (nerror > 0) consensus = best_read_from_graph;
            }
        }
        if (verbosity >=1) {
            cout << "Delete edges between " << left_edge << "  " << right_edge << ": " ;
            edges.Println(cout);
        }
        shb.DeleteEdges(edges);
        int matching_left = 0;
        for(int x = 0; x < min(ltrim, consensus.isize()); x++)
            if (left_seq[x] == consensus[x]) 
                matching_left++;
        int matching_right = 0;
        for(int x = 0; x < min(rtrim, consensus.isize()); x++)
            if (right_seq[rtrim - 1 - x] == consensus[consensus.size() - 1- x]) 
                matching_right++;
        BaseVec new_edge(consensus.begin() + matching_left - K + 1, 
                consensus.end() - matching_right + K - 1); // K-1 overlap with neighbor edges
        if (verbosity >=1) {
            cout << "matching_left= " << matching_left << " ltrim= " << ltrim
                << " matching_right= " << matching_right << " rtrim= " << rtrim 
                << " edge_len= " << new_edge.size() << endl;
        }
        shb.EdgeObjectMutable(left_edge) = BaseVec(ll.begin(), ll.end() - ltrim + matching_left);
        shb.EdgeObjectMutable(right_edge) = BaseVec(rr.begin() + rtrim - matching_right, rr.end());
        // Break the left_edge -> right_edge connection if
        // 1. Ends has been modified
        // Or. New edge different from the (K-1) overlap is added
        int v = to_right[left_edge];
        int w = to_left[right_edge];
        if (v == w) {
            if (rtrim == matching_right && ltrim == matching_left 
                    && consensus == basevector(rr, 0, K-1))
                continue; // don't do anything
            else {
                cout << "*** breaking connection from " << left_edge 
                    << " to " << right_edge << "." << endl;
                // n1 --left--> v(w) --right-> n2
                //    To
                // n1 --left--> v --new edge--> new_node ---right--> n2
                int rightright = to_right[right_edge];
                ForceAssertEq(shb.FromSize(v), 1);
                ForceAssertEq(shb.ToSize(v), 1);
                //ForceAssertEq(shb.ToSize(rightright),1);
                shb.DeleteEdgeFrom(v, 0);

                int new_node = shb.N();
                shb.AddVertices(1);
                shb.AddEdge(v, new_node, new_edge);
                shb.FromMutable(new_node) = {rightright};
                shb.FromEdgeObjMutable(new_node) = {right_edge};
                for(int j = 0; j < shb.ToMutable(rightright).isize(); j++) {
                    if (shb.ToEdgeObjMutable(rightright)[j] == right_edge) {
                        shb.ToMutable(rightright)[j] = new_node;
                        shb.ToEdgeObjMutable(rightright)[j] = right_edge;
                        SortSync(shb.ToMutable(rightright), shb.ToEdgeObjMutable(rightright));
                    }
                }
                status_log << "Deleted " << edges.size() << " edges, new edge added, size= " << 
                    new_edge.size() << endl;
                //DumpGraph(ToString(left_edge) + "_" + ToString(right_edge) + "_added.dot",
                //        shb);
            }
        } else {
            shb.AddEdge(v, w, new_edge);
            status_log << "Deleted " << edges.size() << " edges, new edge added, size= " << 
                new_edge.size();
        }
    }
    // Join continuous edges
    // v1 -> v2 -> v3 -> v4 ====> v1 -> v4
    int njoin = JoinConsecutiveEdges(shb);
    if (verbosity >= 1) cout << "njoin= " << njoin << endl;
    // I don't know why, but sometimes RemoveDeadEdgeObjects crash without disabling inv.
    shb.InvMutable() = vec<int>(shb.EdgeObjectCount(),-1);
    shb.RemoveDeadEdgeObjects( );
    shb.RemoveEdgelessVertices( );
    // do not assign inv
    shb.InvMutable() = vec<int>(shb.EdgeObjectCount(),-1);
    cout << "UnrollWithPacbioRead Summary: " << endl;
    cout << "Total unrolling attemps: " << edge_pairs.size() << endl;
    for (size_t i = 0; i < status_logs.size(); ++i) 
        cout << "Between edge " << edge_pairs[i].first << " and " << edge_pairs[i].second
            << ": " << status_logs[i].str() << endl;
    // Patching gaps
    njoin += PatchWithPacbioReads(longreads, &shb);

    if (njoin >0 && shb.EdgeObjectCount() == 1)
        cout << "Perfect assembly generated after patching using pacbio reads." << endl;

    if (verbosity >= 2) 
        DumpGraph("after3.dot", shb);
}
