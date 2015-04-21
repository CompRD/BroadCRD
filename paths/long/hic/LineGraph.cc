///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Aug 8, 2014 - <crdhelp@broadinstitute.org>
//

#include "paths/long/hic/LineGraph.h"
#include "paths/long/hic/CanonicalTestOrder.h"



void dumpLineGraph( LineGraph const& graph, String const& filename,
        int cseed /* = -1 */, Bool singletons /* = False */,
        Bool edgeLabels  /* = True */)
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
            if ( edgeLabels )
                edge_labels[v].push_back( ToString(obj.count) + " (" +
                        ToString(obj.ratio) + "|" +
                        CanonicalTestOrder::pairIndexToConfigString(obj.alt_config) + ")" );
        }
    }

    if ( !edgeLabels ) edge_labels.clear();

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


#include "graph/DigraphTemplate.h"
template digraphVE<TaggedLineId,EdgeEvidence>::digraphVE();
template int digraphVE<TaggedLineId,EdgeEvidence>::N() const;
template void digraphVE<TaggedLineId,EdgeEvidence>::AddVertex(TaggedLineId const&);
template TaggedLineId const& digraphVE<TaggedLineId,EdgeEvidence>::Vert(int) const;
template void digraphE<EdgeEvidence>::DeleteEdgesFrom( int, const vec<int>& );
template void digraphE<EdgeEvidence>::DeleteEdgesTo( int, const vec<int>& );
template int digraphE<EdgeEvidence>::AddEdge(int, int, EdgeEvidence const&);
template void digraphE<EdgeEvidence>::DeleteEdgesAtVertex(int);
template void digraphVE<TaggedInt,EdgeEvidence>::writeBinary(BinaryWriter&) const;
template void digraphVE<TaggedInt,EdgeEvidence>::readBinary(BinaryReader&);
template Bool digraphE<EdgeEvidence>::TestValid(Bool const) const;
template void digraphE<EdgeEvidence>::Clear();
