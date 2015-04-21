/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "assembly/OffsetGraph.h"

#include "STLExtensions.h"
#include <set>
#include <functional>
#include <iterator>

template <class Entity>
class OffsetGraph<Entity>::Node
{
 public:
  Node( const Entity &entity )
    : m_entity( entity ),
      m_isCollapsed( false )
  {}
  
  void  AddOutEdge( Edge * pEdge, Node *pDestNode )
  {
    m_edgePtrMap.insert( make_pair( pDestNode, pEdge ) );
    m_isCollapsed = false;
  }
  
  Entity  GetEntity() const { return m_entity; }
  
  typedef map<Node*,Edge*> EdgeMap;
  
  typename EdgeMap::const_iterator  BeginEdge() const { return m_edgePtrMap.begin(); }
  typename EdgeMap::const_iterator  EndEdge()   const { return m_edgePtrMap.end(); }
  
  const EdgeMap & GetEdgeMap() const { return m_edgePtrMap; }
  
  Edge* GetBestEdge( OffsetGraph<Entity> &graph, Node* pTarget );
  
 private:
  Entity m_entity;
  
  EdgeMap m_edgePtrMap;
  bool m_isCollapsed;
};

template <class Entity>
class OffsetGraph<Entity>::Edge
{
 public:
  
  Edge( const Offset<Entity> &offset,
        map<Entity,Node*> &nodePtrMap ) 
    : m_offset( offset ),
      m_isBest( false )
  {
    typename map<Entity,Node*>::const_iterator nodeMapIter;
    
    nodeMapIter = nodePtrMap.find( offset.GetAnchored() );
    ForceAssert( nodeMapIter != nodePtrMap.end() );
    m_pOutNode = nodeMapIter->second;
    
    nodeMapIter = nodePtrMap.find( offset.GetFloating() );
    ForceAssert( nodeMapIter != nodePtrMap.end() );
    m_pInNode = nodeMapIter->second;
  }
  
  Offset<Entity>  GetOffset() const { return m_offset; }
  
  void  SetOffset( const Offset<Entity> &newOffset ) 
  {
    Assert( newOffset.GetAnchored() == m_offset.GetAnchored() );
    Assert( newOffset.GetFloating() == m_offset.GetFloating() );
    m_offset = newOffset;
  }
  
  bool  IsBest() const { return m_isBest; }
  
  void  SetIsBest( const bool isBest ) { m_isBest = isBest; }
  
  Node *  GetOutNode() const { return m_pOutNode; }
  Node *  GetInNode()  const { return m_pInNode; }
  
  struct OrderPtrsByDevWithFavored : public binary_function<const Edge*,const Edge*,bool> {
    OrderPtrsByDevWithFavored( Node* pTarget ) 
      : m_pTarget( pTarget ) {}
    bool operator() ( const Edge* pLhs, const Edge* pRhs ) const {
      return ( ( pLhs->GetInNode() != m_pTarget && pRhs->GetInNode() == m_pTarget ) ||
               ( pLhs->GetOffset().GetDev() > pRhs->GetOffset().GetDev() ) );
    }
   private:
    Node* m_pTarget;
  };
  
 private:
  Offset<Entity> m_offset;
  Node *m_pOutNode;
  Node *m_pInNode;
  bool m_isBest;
};


template <class Entity>
void
OffsetGraph<Entity>::AddEdge( const Offset<Entity> &offset )
{
  if ( m_nodePtrMap.count( offset.GetAnchored() ) == 0 )
    m_nodePtrMap.insert( make_pair( offset.GetAnchored(), new Node( offset.GetAnchored() ) ) );
  if ( m_nodePtrMap.count( offset.GetFloating() ) == 0 )
    m_nodePtrMap.insert( make_pair( offset.GetFloating(), new Node( offset.GetFloating() ) ) );
  
  m_edgePtrs.push_back( new Edge( offset, m_nodePtrMap ) );
  m_edgePtrs.back()->GetOutNode()->AddOutEdge( m_edgePtrs.back(),
                                               m_edgePtrs.back()->GetInNode() );
  
  m_edgePtrs.push_back( new Edge( offset.GetSwapped(), m_nodePtrMap ) );
  m_edgePtrs.back()->GetOutNode()->AddOutEdge( m_edgePtrs.back(),
                                               m_edgePtrs.back()->GetInNode() );
}


template <class Entity>
OffsetGraph<Entity>::~OffsetGraph() 
{
  for ( typename map<Entity,Node*>::iterator nodeMapIter = m_nodePtrMap.begin();
        nodeMapIter != m_nodePtrMap.end(); ++nodeMapIter )
    delete nodeMapIter->second;

  for ( typename vec<Edge*>::iterator edgePtrIter = m_edgePtrs.begin();
        edgePtrIter != m_edgePtrs.end(); ++edgePtrIter )
    delete *edgePtrIter;
}


template <class Entity>
void
OffsetGraph<Entity>::AddSuper( const Super &theSuper )
{
  const float numDevs = 3.0;

  WeightedOffsetBuilder<Entity> builder( numDevs );

  typename WeightedOffset<Entity>::Filter theFilter( 0, 0, 
                                                     m_minLinks,
                                                     0, m_minLinks,
                                                     0, m_minLinks,
                                                     INT_MAX );
  builder.SetFilter( theFilter );

  set<Contig> contigsToGraph;
  vec<ContigLocation> contigLocs;
  theSuper.GetContigLocations( contigLocs );
  transform( contigLocs.begin(), contigLocs.end(),
             inserter( contigsToGraph, contigsToGraph.begin() ),
             mem_fun_ref( &ContigLocation::GetContig ) );

  int newSize = m_connected.Size();

  for ( set<Contig>::iterator contigIter = contigsToGraph.begin();
        contigIter != contigsToGraph.end(); ++contigIter )
    if ( m_connectedIdMap.count( *contigIter ) == 0 )
      m_connectedIdMap.insert( make_pair( *contigIter, newSize++ ) );

  if ( m_connected.Size() < newSize ) {
    equiv_rel newConnected( newSize );
    for ( int i = 0; i < m_connected.Size(); ++i )
      newConnected.Join( i, m_connected.ClassId(i) );
    m_connected = newConnected;
  }
  
  vec< WeightedOffset<Entity> > offsets;
  builder.Build( contigsToGraph, offsets );

  for ( unsigned int offsetIdx = 0; offsetIdx < offsets.size(); ++offsetIdx )
  {
    WeightedOffset<Entity> &weightedOffset( offsets[ offsetIdx ] );

    Offset<Entity> offset;
    offset.SetAnchored( weightedOffset.GetAnchored() );
    offset.SetFloating( weightedOffset.GetFloating() );
    offset.SetOffsetAmount( weightedOffset.GetOffsetAmount() );
    offset.SetDev( weightedOffset.GetStdev() );
    offset.SetOrientation( weightedOffset.GetOrient() );

    this->AddEdge( offset );
    
    typename map<Entity,int>::iterator anchoredIdIter = 
      m_connectedIdMap.find( offset.GetAnchored() );
    typename map<Entity,int>::iterator floatingIdIter = 
      m_connectedIdMap.find( offset.GetFloating() );
    
    ForceAssert( anchoredIdIter != m_connectedIdMap.end() );
    ForceAssert( floatingIdIter != m_connectedIdMap.end() );
    m_connected.Join( anchoredIdIter->second, floatingIdIter->second );
  }
}


template <class Entity>
typename OffsetGraph<Entity>::Edge*
OffsetGraph<Entity>::Node::GetBestEdge( OffsetGraph<Entity> &graph, Node* pTarget )
{
  Node *pStartNode = this;

  int maxDev = graph.GetMaxDev();

  typename EdgeMap::const_iterator existingEdgeIter = 
    pStartNode->GetEdgeMap().find( pTarget );

  // cout << "Finding best edge from " 
  //      << pStartNode->GetEntity() << " to " << pTarget->GetEntity() << endl;
  
  if ( existingEdgeIter != pStartNode->GetEdgeMap().end() ) {
    if ( existingEdgeIter->second->IsBest() )
      return existingEdgeIter->second;

    // cout << " Starting with existing edge: " 
    //      << pStartNode->GetEntity() << "\t" 
    //      << existingEdgeIter->second->GetOffset() << "\t" 
    //      << pTarget->GetEntity() << endl;
    
    maxDev = min( maxDev, existingEdgeIter->second->GetOffset().GetDev() );
  }
  
  vec<Edge*> edgesToExamine;
  
  for ( typename EdgeMap::const_iterator firstEdgeIter = pStartNode->BeginEdge();
        firstEdgeIter != pStartNode->EndEdge(); ++firstEdgeIter )
    edgesToExamine.push_back( firstEdgeIter->second );

  //int numEdgesExamined = 0;

  make_heap( edgesToExamine.begin(), edgesToExamine.end(),
             typename Edge::OrderPtrsByDevWithFavored( pTarget ) );
  
  while ( ! edgesToExamine.empty() )
  {
    //Dot( cout, numEdgesExamined++ );

    Edge *pFirstEdge = edgesToExamine.front();

    if ( pFirstEdge->GetOffset().GetDev() > maxDev )
      break;

    pop_heap( edgesToExamine.begin(), edgesToExamine.end(),
              typename Edge::OrderPtrsByDevWithFavored( pTarget ) );
    edgesToExamine.pop_back();
    
    Node *pMidNode = pFirstEdge->GetInNode();
    
    Offset<Entity> firstOffset = pFirstEdge->GetOffset();
    
    for ( typename EdgeMap::const_iterator secondEdgeIter = pMidNode->BeginEdge();
          secondEdgeIter != pMidNode->EndEdge(); ++secondEdgeIter )
    {
      Node *pEndNode = secondEdgeIter->first;
      Edge *pSecondEdge = secondEdgeIter->second;
      
      if ( pEndNode == pStartNode )
        continue;
      
      Offset<Entity> secondOffset = pSecondEdge->GetOffset();
      
      Offset<Entity> jointOffset;
      
      float firstVariance = static_cast<float>( firstOffset.GetDev() );
      firstVariance *= firstVariance;
      
      float secondVariance = static_cast<float>( secondOffset.GetDev() );
      secondVariance *= secondVariance;
      
      int jointDev = static_cast<int>( rintf( sqrtf( firstVariance + secondVariance ) ) );
      
      if ( maxDev >= 0 && jointDev > maxDev )
        continue;
      
      jointOffset.SetAnchored( pStartNode->GetEntity() );
      jointOffset.SetFloating( pEndNode->GetEntity() );
      
      if ( firstOffset.IsForward() )
      {
        jointOffset.SetOffsetAmount( firstOffset.GetOffsetAmount() + 
                                     secondOffset.GetOffsetAmount() );
        jointOffset.SetOrientation( secondOffset.GetOrientation() );
      }
      else
      {
        jointOffset.SetOffsetAmount( firstOffset.GetOffsetAmount() +
                                     pMidNode->GetEntity().GetLength() -
                                     firstOffset.GetOffsetAmount() -
                                     pEndNode->GetEntity().GetLength() );
        jointOffset.SetOrientation( Flip( secondOffset.GetOrientation() ) );
      }
      jointOffset.SetDev( jointDev );

      typename EdgeMap::const_iterator existingEdgeIter = 
        pStartNode->GetEdgeMap().find( pEndNode );
      
      if ( existingEdgeIter != pStartNode->GetEdgeMap().end() )
      {
        if ( jointDev < existingEdgeIter->second->GetOffset().GetDev() )
        {
          // cout << " Replaced existing edge with: " 
          //      << pStartNode->GetEntity() << "\t" 
          //      << jointOffset << "\t" 
          //      << pEndNode->GetEntity() << endl;

          existingEdgeIter->second->SetOffset( jointOffset );
          edgesToExamine.push_back( existingEdgeIter->second );
          push_heap( edgesToExamine.begin(), edgesToExamine.end(),
                     typename Edge::OrderPtrsByDevWithFavored( pTarget ) );
          
          typename EdgeMap::const_iterator reverseEdgeIter;
          reverseEdgeIter = pEndNode->GetEdgeMap().find( pStartNode );
          reverseEdgeIter->second->SetOffset( jointOffset.GetSwapped() );

          if ( pEndNode == pTarget )
            maxDev = min( maxDev, jointOffset.GetDev() );
        }
      }
      else
      {
        // cout << " Added new edge: "
        //      << pStartNode->GetEntity() << "\t" 
        //      << jointOffset << "\t" 
        //      << pEndNode->GetEntity() << endl;
        
        graph.AddEdge( jointOffset );
        
        typename EdgeMap::const_iterator newEdgeIter;
        newEdgeIter = pStartNode->GetEdgeMap().find( pEndNode );
        edgesToExamine.push_back( newEdgeIter->second );
        push_heap( edgesToExamine.begin(), edgesToExamine.end(),
                   typename Edge::OrderPtrsByDevWithFavored( pTarget ) );
        
        if ( pEndNode == pTarget )
          maxDev = min( maxDev, jointOffset.GetDev() );
      }
    }
  }

  // cout << endl;

  // If an edge exists at this point, it's the best one we could find.
  existingEdgeIter = pStartNode->GetEdgeMap().find( pTarget );

  if ( existingEdgeIter != pStartNode->GetEdgeMap().end() ) {
    existingEdgeIter->second->SetIsBest( true );
    
    return existingEdgeIter->second;
  }

  return 0;
}

// Public interface of search method.
template <class Entity>
bool
OffsetGraph<Entity>::FindOffset( const Entity &fromEntity, const Entity &toEntity,
                                 Offset<Entity> &offset )
{
  if ( fromEntity == toEntity ) {
    offset.SetOffsetAmount( 0 );
    offset.SetDev( 0 );
    offset.SetOrientation( orient_FW );
    return true;
  }

  typename map<Entity,int>::const_iterator fromIdIter = m_connectedIdMap.find( fromEntity );
  typename map<Entity,int>::const_iterator toIdIter = m_connectedIdMap.find( toEntity );
  
  ForceAssert( fromIdIter != m_connectedIdMap.end() );
  ForceAssert( toIdIter != m_connectedIdMap.end() );

  if ( ! m_connected.Equiv( fromIdIter->second, toIdIter->second ) ) {
    offset.SetOffsetAmount( 0 );
    offset.SetDev( -1 );
    offset.SetOrientation( orient_FW );
    return false;
  }

  // cout << "Finding path from " << fromEntity << " to " << toEntity << "." << endl;

  typename map<Entity,Node*>::iterator fromIter = m_nodePtrMap.find( fromEntity );
  if ( fromIter == m_nodePtrMap.end() )
    return false;
  Node *pFrom = fromIter->second;

  typename map<Entity,Node*>::iterator toIter = m_nodePtrMap.find( toEntity );
  if ( toIter == m_nodePtrMap.end() )
    return false;
  Node *pTo = toIter->second;

  Edge* pBestEdge = pFrom->GetBestEdge( *this, pTo );

  if ( pBestEdge == 0 )
    return false;

  offset = pBestEdge->GetOffset();

  // cout << "Returning " << offset << endl;

  return true;
}


template OffsetGraph<Contig>::~OffsetGraph();
template void OffsetGraph<Contig>::AddSuper( const Super& );
template bool OffsetGraph<Contig>::FindOffset( const Contig&, const Contig&, Offset<Contig>& );
