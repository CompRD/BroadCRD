/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


#ifndef __INCLUDE_paths_NormalEdgeGraph_h
#define __INCLUDE_paths_NormalEdgeGraph_h

#include <utility>
#include "graph/Digraph.h"
#include "SemanticTypes.h"
#include "system/System.h"
#include "paths/HyperKmerPath.h"
#include "feudal/BinaryStream.h"

// Term: normal edge
//
// A long edge of a HyperKmerPath.   The edges of HyperKmerPaths tend to be
// either very long or very short.  Most of the assembly comes in
// long edges, with tangles of short ones in the middle.


// Semantic type: normal_edge_id_t
// The id of a normal edge; its index in NormalEdgeGraph::normalEdges.
SemanticTypeStd( int, normal_edge_id_t );


/**
   Class: NormalEdgeGraph

   The graph of <normal edges> of a HyperKmerPath.  Each node of the graph
   represents one normal edge.   The edge (n1,n2) is present iff from
   the normal edge n1 you can get to the normal edge n2 without going
   through any other normal edges.

   Often NormalEdgeGraph is passed around not for the links it represents
   between normal edges but simply for the list of normal edges -- the
   correspondence between the normal edge ids and the edge ids in the
   original graph ( OrigEdgeId() and NormalEdgeId() ).
*/
class NormalEdgeGraph: public digraph {

 public:

  NormalEdgeGraph() { }
  NormalEdgeGraph( const HyperKmerPath& hkp, nkmers_t MIN_NORMAL_EDGE_SIZE );

  nkmers_t MinNormalEdgeSize() const { return minNormalEdgeSize; }

  int OrigEdgeId( normal_edge_id_t ne ) const { return normalEdges[ ne ]; }
  normal_edge_id_t NormalEdgeId( int e ) const;

  const vec< int > GetNormalEdges() const { return normalEdges; }

  Bool empty() const { return normalEdges.empty(); }
  Bool nonempty() const { return !empty(); }

  void writeBinary( BinaryWriter& ) const;
  void readBinary( BinaryReader& );
  static size_t externalSizeof() { return 0; }

 private:
  // Field: normalEdges
  // For each normal edge, the edge id of the original edge
  // (in the original HKP).  The normal_edge_id_t identifier of a normal
  // edge gives the normal edge's position in this list.
  vec< int > normalEdges;

  vec< normal_edge_id_t > edge2normal;

  nkmers_t minNormalEdgeSize;

  
};  // class NormalEdgeGraph
SELF_SERIALIZABLE(NormalEdgeGraph);

#endif
// #ifndef __INCLUDE_paths_NormalEdgeGraph_h
