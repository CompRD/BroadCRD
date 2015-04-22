/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


#ifndef ASSEMBLY_OFFSETGRAPH_H
#define ASSEMBLY_OFFSETGRAPH_H

#include "assembly/Offset.h"
#include "Equiv.h"

#include <map>

template <class Entity>
class OffsetGraph
{
 public:
  OffsetGraph()
    : m_minLinks(2),
      m_maxDev(10000)
  {}

  ~OffsetGraph();

  // Contigs connected by fewer links than this are considered unconnected.
  void SetMinLinks( const int minLinks ) { m_minLinks = minLinks; }
  
  int  GetMaxDev( ) const { return m_maxDev; }

  // Offsets with deviation greater than this will not be found.
  void SetMaxDev( const int maxDev ) { m_maxDev = maxDev; }

  // Do not cap offset deviations at all.
  void RemoveDevLimit() { m_maxDev = -1; }

  void  AddSuper( const Super &theSuper );

  // FindOffset(): Finds lowest deviation path from fromEntity to
  // toEntity with deviation less than or equal to maxDev and absolute
  // offset less than or equal to maxAbsOffset, fills out offset of
  // toEntity relative to fromEntity and returns true; if no such path
  // can be found, offset is undefined and returns false.
  bool  FindOffset( const Entity &fromEntity, const Entity &toEntity,
                    Offset<Entity> &offset );

 private:
  class Edge;
  class Node; 

  void  AddEdge( const Offset<Entity> &offset );

  map<Entity,Node*> m_nodePtrMap;
  vec<Edge*>        m_edgePtrs;

  int m_minLinks;
  int m_maxDev;

  equiv_rel m_connected;
  map<Entity,int> m_connectedIdMap;

#if __GNUC__ < 3
  friend class OffsetGraph::Node;
  friend class OffsetGraph::Edge;
#endif
};

typedef OffsetGraph<Contig> ContigOffsetGraph;

#endif
