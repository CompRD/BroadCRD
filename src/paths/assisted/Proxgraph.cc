/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "Superb.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "paths/assisted/CProx.h"
#include "paths/assisted/Proxgraph.h"

/**
 * IsSymmetric
 */
bool IsSymmetric( const proxgraph &pgraph )
{
  // Core info.
  int n_vertices = pgraph.N( );
  int n_edges = pgraph.EdgeObjectCount( );

  vec<int> to_left;
  vec<int> to_right;
  pgraph.ToLeft( to_left );
  pgraph.ToRight( to_right );

  // Create a map (v, w) to (e), where edge e joins v -> w.
  vec< pair<int,int> > vpairs;
  vec<int> eids;
  vpairs.reserve( n_edges );
  eids.reserve( n_edges );
  for (int ii=0; ii<n_edges; ii++) {
    vpairs.push_back( make_pair( to_left[ii], to_right[ii] ) );
    eids.push_back( ii );
  }
  SortSync( vpairs, eids );
  
  // Loop over all connected pairs.
  int n_bad = 0;
  vec< pair<int,int> >::const_iterator it;
  vec< pair<int,int> >::const_iterator vpairs_begin = vpairs.begin( );
  vec< pair<int,int> >::const_iterator vpairs_end = vpairs.end( );

  for (size_t ii=0; ii<vpairs.size( ); ii++) {
    int v = vpairs[ii].first;
    int w = vpairs[ii].second;
    int vrc = ( v % 2 == 0 ) ? v + 1 : v - 1;
    int wrc = ( w % 2 == 0 ) ? w + 1 : w - 1;
    pair<int,int> rcpair = make_pair( wrc, vrc );
    it = lower_bound( vpairs_begin, vpairs_end, rcpair );
    ForceAssert( ! ( it == vpairs_end ) );
    
    if ( ! ( *it == rcpair ) ) n_bad++;
  }
  
  // Done.
  return ( n_bad == 0 );
}

/**
 * AreIdentical
 */
bool AreIdentical( const proxgraph &left,
		   const proxgraph &right,
		   ostream *log )
{
  ofstream devnull ( "/dev/null" );
  ostream &out = log ? *log : devnull;

  out << "Starting to test proxgraphs identity" << endl;

  out << "Testing vertices count: ";
  if ( left.N( ) != right.N( ) ) {
    out << "FAILED (left has " << left.N( )
	<< " vertices, right " << right.N( )
	<< ")" << endl;
    return false;
  }
  out << "OK" << endl;
  
  out << "Testing edge count: ";
  if ( left.EdgeObjectCount( ) != right.EdgeObjectCount( ) ) {
    out << "FAILED (left has " << left.EdgeObjectCount( )
	<< " edges, right " << right.EdgeObjectCount( )
	<< ")" << endl;
    return false;
  }
  out << "OK" << endl;

  out << "Testing " << ToStringAddCommas(left.EdgeObjectCount()) << " edges: ";
  vec<int> to_leftL;
  vec<int> to_rightL;
  left.ToLeft( to_leftL );
  left.ToRight( to_rightL );
  vec<int> to_leftR;
  vec<int> to_rightR;
  right.ToLeft( to_leftR );
  right.ToRight( to_rightR );
  for (int ii=0; ii<left.EdgeObjectCount( ); ii++) {
    int vL = to_leftL[ii];
    int wL = to_rightL[ii];
    int vR = to_leftR[ii];
    int wR = to_rightR[ii];
    if ( vL != vR || wL != wR ) {
      out << "FAILED (edge " << ii
	  << " goes " << vL
	  << " -> " << wL
	  << " in left, but " << wR
	  << " -> " << wR
	  << " in right)" << endl;
      return false;
    }
  }
  out << "OK" << endl;

  out << "Proxgraph identity test passed" << endl;
  return true;
}

/**
 * RemoveContraryRefLinks: remove edges from verticies which are
 * supported by reference only where there are other edges supported
 * by links.
 */
void RemoveContraryRefLinks(proxgraph &pgraph)
{
  int N = pgraph.N();
  for (int v = 0; v < N; ++v) {
    vec<int> links = pgraph.ToEdgeObj(v);
    bool jumpLinks = false;
    vec<int> refonly;
    vec<int> jumponly;
    for (int i = 0; i < links.isize(); ++i) {
      int link = links[i];
      CProx prox = pgraph.EdgeObject(link);
      // keep track of unconfirmed reference or jump links
      if (prox.IsRefGap() && !prox.IsLinkGap())
	refonly.push_back(link);
      if (prox.IsLinkGap() && !prox.IsRefGap())
	jumponly.push_back(link);
    }
    if (refonly.size() > 0 && jumponly.size() > 0) {
      for (int i = 0; i < refonly.isize(); ++i) {
	cout << "deleteing contrary ref edge " << refonly[i] << " from " << v << endl; 
      }
      pgraph.DeleteEdges(refonly);
    }
  }
}

/**
 * Template instantiations
 */
#include "graph/DigraphTemplate.h"

template proxgraph::digraphE();
template proxgraph::digraphE( const vec< vec<int> >&,
			      const vec< vec<int> >&,
			      const vec<CProx>&,
			      const vec< vec<int> >&,
			      const vec< vec<int> >&, const Bool );
template void proxgraph::DumpGraphML( const String& ) const;
template const CProx& proxgraph::EdgeObject( int ) const;
template int proxgraph::EdgeObjectCount() const;
template int proxgraph::EdgeObjectIndexByIndexFrom(int, int) const;
template int proxgraph::EdgeObjectIndexByIndexTo(int, int) const;
template vec<int> const& proxgraph::FromEdgeObj(int) const;
template void proxgraph::Initialize(vec<vec<int> > const&,
                                                    vec<vec<int> > const&,
                                                    vec<CProx> const&,
                                                    vec<vec<int> > const&,
                                                    vec<vec<int> > const&, const Bool);
template int proxgraph::AddEdge(int, int, const CProx &);
template vec<int> const& proxgraph::ToEdgeObj( int ) const;
template vec<int> proxgraph::EdgesBetween( const int, const int ) const;
template void proxgraph::SetEdgeObject( int, const CProx&);
template void proxgraph::ToLeft ( vec<int>& ) const;
template void proxgraph::ToRight ( vec<int>& ) const;
template void proxgraph::readBinary( BinaryReader& );
template void proxgraph::writeBinary( BinaryWriter& ) const;
template void proxgraph::DeleteEdges(const vec<int> &);
template void proxgraph::DeleteEdgeFrom(int, int);
template int proxgraph::InputFromOutputTo(int, int) const;
