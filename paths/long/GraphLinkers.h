///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#ifndef GRAPH_LINKERS_H
#define GRAPH_LINKERS_H

#include "Basevector.h"
#include "PairsManager.h"
#include "Qualvector.h"
#include "paths/HyperBasevector.h"
#include "paths/UnibaseUtils.h"
#include "kmers/KMer.h"

//  Find the read pairs that link two graph components. Only lookup the first
//  kmer of the reads on the graphs. 
//
//  The implementation is just for proof-of-principle, and is extremely
//  inefficient.

class GraphLinkers {
public:
    GraphLinkers(const String& head) 
        :reads(head+".fastb"), quals(head+".qualb"), pm(head + ".pairs")
    {};

    void FindLinks(int K, const HyperBasevector& hbv, const String& dot);

    void FilterLinks();

    // Print all the links between compoments.
    void PrintLinks(const HyperBasevector& hbv) const;

    // Print all reads linking two edges
    void ShowReads(int edge1, int edge2, const HyperBasevector& hbv);

    // Show links between the hyperbasevector edges in dot file. 
    void VisualizeLinks(const HyperBasevector& hbv, const String& dot);

    // Return a graph patched using the links. 
    // !!! Note that current implementation create zero-length edges for links.
    HyperBasevector PatchedGraph(const HyperBasevector& hbv) const;

    // Visualize a patched graph
    static void VisualizePatched(const HyperBasevector& hbv, const String& dot);

    static void SimplifyPatched(HyperBasevector& hbv);

private:
    typedef tuple<int,int,int,int,int,int> Link; 
    typedef vec<Link> LinksVec; 

    // Internal function that will dispatch the call FindLinks using different Ks.
    template <int K>
    void FindLinksK(const HyperBasevector& hbv);

    // A link between two edges is represented as  a tuple of { edge1, edge2,
    // rid1, rid2, d1, d2 }
    //
    //          |-- d1 ---|     |-- d2 --|
    //      rid1--->                 <---- rid2
    //      ---(edge1)----       --------(edge2)-------
    //
    LinksVec links; 
    vecbvec reads;
    vecqvec quals;
    PairsManager pm;
};

#endif

