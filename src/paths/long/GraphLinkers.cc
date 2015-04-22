///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
// MakeDepend: cflags OMP_FLAGS
#include "paths/long/GraphLinkers.h"
#include "PairsManager.h"
#include "paths/HyperBasevector.h"
#include <cctype>

namespace {

// Assign component id to each each of the graph
vec<int> ComponentID(const HyperBasevector& hbv)
{
    vec<int> component_ids(hbv.EdgeObjectCount(),-1);
    vec<vec<int>> components;
    hbv.ComponentEdges(components);
    for (size_t i = 0; i < components.size(); i++) 
        for (int eid: components[i]) 
            component_ids[eid] = i;
    return component_ids;
}

}

// Internal function that will dispatch the call FindLinks using different Ks.
template <int K>
void GraphLinkers::FindLinksK(const HyperBasevector& hbv)
{
    cout << Date() << ": Finding graph linkers. " << endl;
    bool BoundaryEdgeOnly = true;
    bool ExtendEdge       = true;    // Extend each edgep 1000 base upstream
    int  ExtendLen        = 500;

    cout << Date() << ": Building kmer lookup table." << endl;
    typedef KMer<K> Kmer;
    vec<int> component_ids = ComponentID(hbv);
    vec<int> to_right; hbv.ToRight(to_right);
    vec<int> to_left; hbv.ToLeft(to_left);
    // find all edges to build kmer lookups
    vec<int> end_edges;
    for (int i = 0; i < hbv.EdgeObjectCount(); i++) {
        int vert_right = to_right[i];
        if (hbv.FromSize(vert_right) == 0 || !BoundaryEdgeOnly) 
            end_edges.push_back(i);
    }
    cout << "end_edges.size()= " << end_edges.size() << endl;
    // kmer lookup in the edges
    vec<tuple<Kmer, int, int>> kmer_edge_offset;
    #pragma omp parallel
    {
        vec<tuple<Kmer, int, int>> kplus_local;
        #pragma omp for
        for (size_t i = 0; i < end_edges.size(); i++) {
            int ee = end_edges[i];
            int vert_right = to_right[ee];

            vec<bool> visited(hbv.EdgeObjectCount(), false);
            vec<tuple<int,int, int>> todo; // edge_id, kmer_start, dist_to_end
            todo.push(ee,  hbv.EdgeLengthBases(ee) - K, K);
            visited[ee] = true;
            int edge, pos, to_end;
            while (!todo.empty()) {
                tie(edge, pos, to_end) = todo.back();
                todo.pop_back();
                if (to_end > ExtendLen) continue;
                const basevector& base = hbv.EdgeObject(edge);
                while (pos >= 0 && to_end < ExtendLen) {
                    kplus_local.push(Kmer(base.begin()+pos), ee, to_end);
                    --pos,  ++to_end;
                }
                if (pos == -1 && to_end < ExtendLen) {
                    int vert_left = to_left[edge];
                    for (int j = 0; j < hbv.ToSize(vert_left); j++) {
                        int edge2 = hbv.EdgeObjectIndexByIndexTo(vert_left, j);
                        if (!visited[edge2]) {
                            todo.push(edge2,  hbv.EdgeLengthBases(edge2) - K, to_end);
                            visited[edge2] = true;
                        }
                    }
                }
            }
        }
        #pragma omp critical
        kmer_edge_offset.append(kplus_local);
    }
    Sort(kmer_edge_offset);
    cout << Date() << ": Done creating kmer lookup table. " << endl;
    // edge and reverse complement
    vecbasevector edges;
    for (int i = 0; i < hbv.EdgeObjectCount(); i++) 
        edges.push_back(hbv.EdgeObject(i));
    vec<int> to_rc;
    UnibaseInvolution(edges, to_rc);

    cout << Date() << ": Find read pair locations on edges." << endl;
    #pragma omp parallel for
    for ( size_t pairId = 0; pairId < pm.nPairs(); ++pairId ){
        auto FindEdgesHomes = [&](int rid) {
            const int MinQual = 3;
            int len = 0;
            while (len != (int)this->quals[rid].size() && this->quals[rid][len] >= MinQual)
                len++;
            vec<pair<int,int>> edges;
            for (int i = 0; i < len - K + 1; i++) {
                Kmer kmer(this->reads[rid].begin() + i);
                auto it1 = lower_bound(kmer_edge_offset.begin(), kmer_edge_offset.end(),
                        make_tuple(kmer, 0, 0));
                while (it1 != kmer_edge_offset.end() && get<0>(*it1) == kmer) {
                    edges.push(get<1>(*it1), get<2>(*it1) + i);
                    ++it1;
                }
                // Only one kmer
                break;
            }
            return std::move(edges);
        };

        int rid1 = pm.ID1(pairId);
        int rid2 = pm.ID2(pairId);
        vec<pair<int,int>> edge1s = FindEdgesHomes(rid1);
        vec<pair<int,int>> edge2s = FindEdgesHomes(rid2);
        //if (edge1s.empty() || edge2s.empty()) continue;

        set<Link> link_pool; // links can be duplcate more than head kmers are used
        for (const pair<int,int>& p1: edge1s)
        for (const pair<int,int>& p2: edge2s) {
            int edge1    = p1.first;
            int d1       = p1.second;
            int edge1_rc = to_rc[edge1];
            int edge2    = p2.first;
            int d2       = p2.second;
            int edge2_rc = to_rc[edge2];
            if (component_ids[edge1] != component_ids[edge2_rc])
                link_pool.insert(Link(edge1, edge2_rc, rid1, rid2, p1.second, p2.second));
            if (component_ids[edge1_rc] != component_ids[edge2]) 
                link_pool.insert(Link(edge2, edge1_rc, rid2, rid1, p2.second, p1.second));
        }
        #pragma omp critical
        copy(link_pool.begin(), link_pool.end(), back_inserter(links));
    }
    Sort(links);
}

void GraphLinkers::FindLinks(int K, const HyperBasevector& hbv, const String&
        dot) 
{
    switch (K) {
        case 21: FindLinksK<21>(hbv); break;
        case 31: FindLinksK<31>(hbv); break;
        case 41: FindLinksK<41>(hbv); break;
        case 51: FindLinksK<51>(hbv); break;
        case 61: FindLinksK<61>(hbv); break;
        default: cout << "Unsupported K " << K << endl;
                 exit(1);
    }
    FilterLinks();
    //PrintLinks(hbv);
    //VisualizeLinks(hbv, dot);
}

// Print lists of links between graph components. 

void GraphLinkers::PrintLinks(const HyperBasevector& hbv) const
{
    vec<int> component_ids = ComponentID(hbv);;
    std::map<pair<int,int>, LinksVec> links_grouped;
    for (auto& x: links) {
        int group1 = component_ids[get<0>(x)];
        int group2 = component_ids[get<1>(x)];
        links_grouped[make_pair(group1, group2)].push_back(x);
    }
    for (auto it = links_grouped.begin(); it != links_grouped.end(); ++it) {
        int nlinks = it->second.size();
        cout << nlinks << " between components " << it->first.first << " and " 
            << it->first.second << endl;
        for (auto x: it->second) {
            int edge1, edge2, rid1, rid2, to_end1, to_end2;
            tie(edge1, edge2, rid1, rid2, to_end1, to_end2) = x;
            cout << "    edges " << edge1 << " " << edge2 << " by pair " 
                << rid1 << " " << rid2  
                << " sep = L - " << to_end1 + to_end2 << endl;;
        }
    }
}

// Create a dot file to display the links between the graph components.
// The links between edges are shown in dotted lines.

void GraphLinkers::VisualizeLinks(const HyperBasevector& hbv, const String& out_head)
{
    const int AvgFragLen = 400;
    vec<tuple<int,int,int,int>> bridge; // edge1, edge2, nlinks, sep
    size_t i = 0, j = 0;
    while (i < links.size()) {
        int e1 = get<0>(links[i]);
        int e2 = get<1>(links[i]);
        int x1 = get<4>(links[i]);
        int x2 = get<5>(links[i]);
        j = i;
        vec<int> seps;
        while (j < links.size() && get<0>(links[j]) == e1 
                && get<1>(links[j]) == e2) {
            seps.push_back(AvgFragLen - get<4>(links[j]) - get<5>(links[j]));
            j++;
        }
        Sort(seps);
        int sep_med = Median(seps);
        int n = j - i;
        bridge.push(e1, e2, n, sep_med);
        i = j;
    }

    HyperBasevector hb = hbv;
    vec<int> to_left, to_right;
    hb.ToLeft(to_left); 
    hb.ToRight(to_right);
    int nedges_total = hb.EdgeObjectCount() + bridge.size();
    vec<double> length(nedges_total, 0);
    vec<String> labels(nedges_total,"");
    vec<Bool> dashed(nedges_total, False);
    for (int i = 0; i < hb.EdgeObjectCount(); i++) {
        length[i] = hb.EdgeLengthBases(i);
        labels[i] = ToString(i);
    }
    for (auto& x: bridge) {
        int e1, e2, n, sep;
        tie(e1, e2, n, sep) = x;
        int left = to_right[e1];
        int right = to_left[e2];
        hb.AddEdge(left, right, basevector("")); // empty edge
        int edge_index = hb.EdgeObjectCount() - 1;
        labels[edge_index] = ToString(n) + "(" + ToString(sep) + ")";
        length[edge_index] = sep;
        dashed[edge_index] = True;
    }
    ofstream out(out_head + ".linked.dot");
    hb.PrettyDOT(out, length, HyperBasevector::edge_label_info(HyperBasevector::
                edge_label_info::DIRECT, &labels), 
            True, False, NULL, NULL, NULL, &dashed);

    //ofstream out2(out_head + ".unlinked.dot");
    //hbv.PrintSummaryDOT0w(out2);
}

void GraphLinkers::FilterLinks()
{
    const int MaxFragLen = 600;
    const int MinLinks = 2;
    vec<Bool> todel(links.size(), false);
    size_t i = 0, j = 0;
    while (i < links.size()) {
        int e1 = get<0>(links[i]);
        int e2 = get<1>(links[i]);
        int x1 = get<4>(links[i]);
        int x2 = get<5>(links[i]);
        j = i;
        int n_valid_links = 0;
        while (j < links.size() && get<0>(links[j]) == e1 
                && get<1>(links[j]) == e2) {
            int d1 = get<4>(links[j]);
            int d2 = get<5>(links[j]);
            if (d1 + d2 < MaxFragLen) 
                n_valid_links++;
            else
                todel[j] = true;
            j++;
        }
        if (n_valid_links < MinLinks)
            fill(todel.begin() + i, todel.begin() + j, true);
        i = j; 
    }
    EraseIf(links, todel);
}


void GraphLinkers::ShowReads(int edge1, int edge2, const HyperBasevector& hbv) 
{
    const int MinQual = 3;
    cout << Date() << ": Show reads between edges " << edge1 << " " << edge2 << endl;
    for (int edge: {edge1, edge2})
        if (edge < 0 || edge >= hbv.EdgeObjectCount()) {
            cout << "Invalid edge id " << edge << endl;
            CRD::exit(1);
        }
    vec<pair<int,int>> reads_on_edge1, reads_on_edge2;
    auto it1 = lower_bound(links.begin(), links.end(), 
            make_tuple(edge1, edge2, -1e9, -1e9, 0, 0));
    auto it2 = upper_bound(links.begin(), links.end(), 
            make_tuple(edge1, edge2, 1e9, 1e9, 0, 0));
    for (auto it = it1; it != it2; ++it) {
        int e1, e2, rid1, rid2, d1, d2;
        tie(e1, e2, rid1, rid2, d1, d2) = *it;
        if (e1 == edge1 && e2 == edge2) {
            reads_on_edge1.push(rid1, -d1);
            reads_on_edge2.push(rid2, d2 - reads[rid2].isize());
        }
    }
    UniqueSort(reads_on_edge1);
    UniqueSort(reads_on_edge2);

    int min_offset = -hbv.EdgeLengthBases(edge1);
    for (auto x: reads_on_edge1) min_offset = min(min_offset, x.second);
    for (auto x: reads_on_edge2) min_offset = min(min_offset, x.second);

    int actual_padding = -min_offset;
    vec<String> lines;
   
    String line =  String(actual_padding - hbv.EdgeLengthBases(edge1), ' ') 
            + hbv.EdgeObject(edge1).ToString() + " edge_" + ToString(edge1);
    lines.push_back(line);
    for (size_t i = 0; i < reads_on_edge1.size(); i++) {
        int rid1 = reads_on_edge1[i].first;
        String read_str = reads[rid1].ToString();
        for (size_t j = 0; j < read_str.size(); j++) {
            const int MinQual = 3;
            if (quals[rid1][j] < MinQual)
                read_str[j] = tolower(read_str[j]);
        }
        line = String(reads_on_edge1[i].second + actual_padding, ' ')
                + read_str + " read_" + ToString(rid1);
        lines.push_back(line);
    }

    String padding(actual_padding, ' ');
    String label = "edge_" + ToString(edge2) + " ";
    if (padding.size() >= label.size())
        copy(label.begin(), label.end(), padding.end() - label.size());
    line =  padding + hbv.EdgeObject(edge2).ToString();
    lines.push_back(line);
    for (size_t i = 0; i < reads_on_edge2.size(); i++) {
        int rid2 = reads_on_edge2[i].first;
        String padding(actual_padding + reads_on_edge2[i].second, ' ');
        String label = "read_" + ToString(rid2) + " ";
        if (padding.size() >= label.size())
            copy(label.begin(), label.end(), padding.end() - label.size());

        String read_str = ReverseComplement(reads[rid2]).ToString();
        for (size_t j = 0; j < read_str.size(); j++)
            if (quals[rid2][read_str.size() - 1 - j] < MinQual)
                read_str[j] = tolower(read_str[j]);
        line = padding + read_str;
        lines.push_back(line);
    }

    const int LeftWindow = 150;
    int trim_left = actual_padding - LeftWindow;
    for (auto &x: lines) trim_left = min(trim_left, (int)x.size());
    for (size_t i = 0; i < lines.size(); i++) {
        if (trim_left > 0) {
            cout << String(lines[i].begin()+trim_left, lines[i].end()) << endl;
        }
        else
            cout << String(-trim_left, ' ') << lines[i] << endl;
    }
}

// Very crude way of patching graphs. The edges are simply linked by zero-length
// basevector. Note that it violates the K-1 overlap assumption for normal
// Hyperbasevector object.

HyperBasevector GraphLinkers::PatchedGraph(const HyperBasevector& hbv) const
{
    int verbosity = 0;

    cout << Date() << ": Before patching, there are " << hbv.NComponents() 
        << " components" << endl;
    vec<pair<int,int>> bridge_pair;
    vec<int> bridge_nlinks;
    for (size_t i = 0; i < links.size(); i++) {
        size_t j = i+1;
        int e0 = get<0>(links[i]), e1 = get<1>(links[i]);
        while (j < links.size() && get<0>(links[j]) == e0 && 
                    get<1>(links[j]) == e1) 
            j++;
        bridge_pair.push(e0, e1);
        bridge_nlinks.push(j-i);
        i = j-1;
    }
    ReverseSortSync(bridge_nlinks, bridge_pair);

    HyperBasevector hb = hbv;
    vec<int> to_left, to_right;
    hb.ToLeft(to_left); 
    hb.ToRight(to_right);

    for (size_t i = 0; i < bridge_pair.size(); i++) {
        int e0 = bridge_pair[i].first;
        int e1 = bridge_pair[i].second;
        int left = to_right[e0];
        int right = to_left[e1];
        if (hb.FromSize(left) > 0 || hb.ToSize(right)) {
            if (verbosity >= 1)
                cout << "Skip link between " << e0 << " " << e1 
                    << " nlinks= " << bridge_nlinks[i] << endl;
        }
        hb.AddEdge(left, right, basevector("")); // empty edge
    }
    cout << Date() << ": After patching, there are " << hb.NComponents() 
        << " components" << endl;
    return hb;
}


void GraphLinkers::VisualizePatched(const HyperBasevector& hbv, const String& dot)
{
    int nedge = hbv.EdgeObjectCount();
    vec<double> length(nedge, 0);
    vec<String> labels(nedge,"");
    vec<Bool> dashed(nedge, False);
    for (int i = 0; i < nedge; i++) {
        length[i] = hbv.EdgeLengthBases(i);
        labels[i] = ToString(i);
        dashed[i] = (length[i] == 0);
    }
    ofstream out(dot);
    hbv.PrettyDOT(out, length, HyperBasevector::edge_label_info(HyperBasevector::
                edge_label_info::DIRECT, &labels), 
            True, False, NULL, NULL, NULL, &dashed);
}


void GraphLinkers::SimplifyPatched(HyperBasevector& hbv) 
{
    // Remove Dangle
    // ----->(v2)------->(v)---
    //         \--(e1)-->(v1)
    const int MinLen = 100;
    for (bool changed = true; changed; ){
        changed = false;
        for (int v1 = 0; v1 < hbv.N(); v1++) {
            if (hbv.FromSize(v1) == 1 && hbv.ToSize(v1) == 0) {
                int e1 = hbv.EdgeObjectIndexByIndexFrom(v1, 0);
                int v2 = hbv.From(v1)[0];
                if (hbv.ToSize(v2) + hbv.FromSize(v2) > 2 &&
                        hbv.EdgeLengthBases(e1) < MinLen) {
                    hbv.DeleteEdgeFrom(v1, 0);
                    changed = true;
                    break;
                }
            }
            if (hbv.FromSize(v1) == 0 && hbv.ToSize(v1) == 1) {
                int e1 = hbv.EdgeObjectIndexByIndexTo(v1, 0);
                int v2 = hbv.To(v1)[0];
                if (hbv.ToSize(v2) + hbv.FromSize(v2) > 2 &&
                        hbv.EdgeLengthBases(e1) < MinLen) {
                    hbv.DeleteEdgeTo(v1, 0);
                    changed = true;
                    break;
                }
            }
        }
    } 

    // Falttern the bubbles
    // -> v1 ----e0---->v2 --  ===>  -- v1 ---e1--> v2 ---
    //     | ----e1 --->|
    //      -----e2 --->
    const int MaxLenthDiff = 100;
    for (int v1 = 0; v1 < hbv.N(); v1++) {
        if (hbv.ToSize(v1) != 1) continue;
        if (hbv.FromSize(v1) < 2) continue;
        int e0 = hbv.EdgeObjectIndexByIndexFrom(v1, 0);
        int v2 = hbv.From(v1)[0];
        int len0 = hbv.EdgeLengthBases(e0);
        bool ends_match = true;
        bool length_match = true;
        bool has_patch = (hbv.EdgeLengthBases(e0) == 0);
        for (int j = 1; j < hbv.FromSize(v1); j++) {
            if (hbv.From(v1)[j] != v2){
                ends_match = false; break;
            }
            int e = hbv.EdgeObjectIndexByIndexFrom(v1, j);
            if (hbv.EdgeLengthBases(e) == 0){
                has_patch = true;
                break;
            }
            if (abs(hbv.EdgeLengthBases(e)-len0) > MaxLenthDiff) {
                length_match = false;
                break;
            }
        }
        if (hbv.ToSize(v2) == hbv.FromSize(v1) && ends_match && length_match
                && !has_patch) {
            for (int j = 1; j < hbv.FromSize(v1); j++) {
                hbv.DeleteEdgeFrom(v1, j);
            }
        }
    }

    // join edges (v1)--e1-->(v2)--e2-->(v3) ===> v1 ---e3---> v3
    for (bool changed = true; changed; ) {
        changed = false;
        for (int v2 = 0; v2 < hbv.N(); v2++) {
            if (hbv.ToSize(v2) != 1) continue;
            if (hbv.FromSize(v2) != 1) continue;
            int v3 = hbv.From(v2)[0];
            int v1 = hbv.To(v2)[0];
            int e12 = hbv.EdgeObjectIndexByIndexTo(v2, 0);
            int e23 = hbv.EdgeObjectIndexByIndexFrom(v2, 0);
            if (e12 == e23) continue;
            if (hbv.EdgeLengthBases(e12) == 0 ||
                    hbv.EdgeLengthBases(e23) == 0) continue;
            vec<int> path = {e12, e23};
            basevector e13 = hbv.EdgePathToBases(path);
            hbv.DeleteEdgeTo(v2, 0);
            hbv.DeleteEdgeFrom(v2, 0);
            hbv.AddEdge(v1, v3, e13);
            changed = true;
            break;
        }
    } 

    hbv.RemoveSmallComponents(100);
}

