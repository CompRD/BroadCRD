///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "paths/assisted/Backbone.h"
#include "paths/assisted/CProx.h"
#include "paths/assisted/Proxgraph.h"

/**
 * BackboneFromCProx
 *
 * It takes in input a proxgraph and it generates in output a
 * trivial backbone digraphVE (see Backbone.h for details).
 */
void BackboneFromCProx( const proxgraph &in_graph,
			backbone &out_graph )
{
  int n_vertices = in_graph.N( );
  vec< vec<int> > vert( n_vertices, vec<int>( 1, 0 ) );
  for (int ii=0; ii<n_vertices; ii++)
    vert[ii][0] = ii;
  
  backbone graph( in_graph, vert );
  swap( graph, out_graph );
}

/**
 * BackboneToCProx
 *
 * It generates a proxgraph from a given backbone. If in_graph
 * is a proxgraph, and if in_backbone is the trivial backbone
 * generated from in_graph, then out_graphE is guaranteed to match
 * exactly in_graphE.
 */
void BackboneToCProx( const backbone &in_backbone,
		      const proxgraph &in_graphE,
		      proxgraph &out_graphE )
{
  int n_vertices = in_graphE.N( );
  int n_backbone_vertices = in_backbone.N( );
  int n_backbone_edges = in_backbone.EdgeObjectCount( );
  
  // Core descriptor for out_graphE.
  vec< vec<int> > from( n_vertices );
  vec< vec<int> > to( n_vertices );
  vec< vec<int> > fromEO( n_vertices );
  vec< vec<int> > toEO( n_vertices );
  vec<CProx> edges;

  // after[v]:
  //  w    if v -> w in some V in in_backbone
  // -1    all other cases
  vec<int> after( n_vertices, -1 );
  for (int ii=0; ii<n_backbone_vertices; ii++) {
    const vec<int> &vert = in_backbone.Vert( ii );
    for (int jj=0; jj<vert.isize( )-1; jj++)
      after[ vert[jj] ] = vert[jj+1];
  }
  
  // Standard to_left, to_right.
  vec<int> to_left;
  vec<int> to_right;
  in_graphE.ToLeft( to_left );
  in_graphE.ToRight( to_right );

  vec<int> to_left_bb;
  vec<int> to_right_bb;
  in_backbone.ToLeft( to_left_bb );
  in_backbone.ToRight( to_right_bb );
  
  // First add edges between consecutive vertices in each V in in_backbone.
  int n_edges_in_graphE = in_graphE.EdgeObjectCount( );
  for (int eid=0; eid<n_edges_in_graphE; eid++) {
    int vid = to_left[eid];
    int wid = to_right[eid];
    if ( wid != after[vid] ) continue;
    
    from[vid].push_back( wid );
    fromEO[vid].push_back( edges.isize( ) );
    to[wid].push_back( vid );
    toEO[wid].push_back( edges.isize( ) );
    edges.push_back( in_graphE.EdgeObject( eid ) );
  }

  // Then edges between V in in_backbone.
  for (int eid=0; eid<n_backbone_edges; eid++) {
    int Vid = to_left_bb[eid];
    int Wid = to_right_bb[eid];
    if ( Vid == Wid ) continue;

    const vec<int> &vids = in_backbone.Vert( Vid );
    const vec<int> &wids = in_backbone.Vert( Wid );
    int vid = vids[vids.size( )-1];
    int wid = wids[0];
    
    from[vid].push_back( wid );
    fromEO[vid].push_back( edges.isize( ) );
    to[wid].push_back( vid );
    toEO[wid].push_back( edges.isize( ) );
    edges.push_back( in_backbone.EdgeObject( eid ) );
  }

  // Sort all.
  for (int ii=0; ii<n_vertices; ii++) {
    sort( from[ii].begin( ), from[ii].end( ) );
    sort( to[ii].begin( ), to[ii].end( ) );
    sort( fromEO[ii].begin( ), fromEO[ii].end( ) );
    sort( toEO[ii].begin( ), toEO[ii].end( ) );
  }

  // Build graph.
  out_graphE.Initialize( from, to, edges, toEO, fromEO );

}


/**
 * IsSymmetric
 */
bool IsSymmetric( const backbone &bbg, vec<int> *VertRC, vec<int> *EdgeRC )
{
  // Local rc maps.
  vec<int> toVRC;
  vec<int> toERC;
  
  // Core info.
  int n_Vertices = bbg.N( );
  int n_edges = bbg.EdgeObjectCount( );
  int n_vertices = 0;
  for (int ii=0; ii<n_Vertices; ii++)
    if ( bbg.Vert( ii ).size( ) > 0 )
      n_vertices = Max( n_vertices, Max( bbg.Vert( ii ) ) );
  n_vertices++;

  vec<int> to_left;
  vec<int> to_right;
  bbg.ToLeft( to_left );
  bbg.ToRight( to_right );

  // Map vertices to Vertices (and check all vertices are accounted for).
  vec<int> toV( n_vertices, -1 );
  for (int ii=0; ii<n_Vertices; ii++) {
    const vec<int> &Vert = bbg.Vert( ii );
    for (int jj=0; jj<Vert.isize( ); jj++)
      toV[ Vert[jj] ] = ii;
  }
  for (int ii=0; ii<n_vertices; ii++)
    ForceAssertNe( toV[ii], -1 );

  // Try to build the map VertRC. Will check it later.
  toVRC.resize( n_Vertices, -1 );
  for (int ii=0; ii<n_Vertices; ii++) {
    const vec<int> &Vert = bbg.Vert( ii );
    if ( Vert.size( ) < 1 ) continue;
    int v = Vert[0];
    int vrc = ( v % 2 == 0 ) ? v + 1 : v - 1;
    int iiRC = toV[vrc];
    toVRC[ii] = iiRC;
    toVRC[iiRC] = ii;
  }

  // Check toVRC.
  bool ok = true;
  for (int ii=0; ii<toVRC.isize( ); ii++) {
    if ( ! ok ) break;

    int jj = toVRC[ii];
    ForceAssertNe( ii, jj );   // V cannot be its own rc
    if ( jj < ii ) continue;   // avoid duplication (and skip empty Vertices)
    
    const vec<int> &Vert = bbg.Vert( ii );
    const vec<int> &VertRC = bbg.Vert( jj );
    if ( Vert.size( ) != VertRC.size( ) ) 
      ok = false;
    else {
      for (int kk=0; kk<Vert.isize( ); kk++) {
	int v = Vert[kk];
	int w = VertRC[Vert.isize( ) - 1 - kk];
	if ( v != w - 1 && v != w + 1 ) {
	  ok = false;
	  break;
	}
      }
    }
  }
  
  // Failed.
  if ( ! ok ) {
    toVRC.clear( );
    return false;
  }
  
  // Create a map (V, W) to (e), where edge e joins V -> W.
  vec< pair<int,int> > Vpairs;
  vec<int> eids;
  Vpairs.reserve( n_edges );
  eids.reserve( n_edges );
  for (int ii=0; ii<n_edges; ii++) {
    ForceAssert( bbg.Vert( to_left[ii] ).size( ) > 0 );
    ForceAssert( bbg.Vert( to_right[ii] ).size( ) > 0 );
    Vpairs.push_back( make_pair( to_left[ii], to_right[ii] ) );
    eids.push_back( ii );
  }
  SortSync( Vpairs, eids );
  
  // Check edges, and fill toERC.
  toERC.resize( n_edges, -1 );
  
  // Loop over all connected pairs.
  int n_bad = 0;
  vec< pair<int,int> >::const_iterator it;
  vec< pair<int,int> >::const_iterator itRC;
  vec< pair<int,int> >::const_iterator Vpairs_begin = Vpairs.begin( );
  vec< pair<int,int> >::const_iterator Vpairs_end = Vpairs.end( );

  for (size_t ii=0; ii<Vpairs.size( ); ii++) {
    int V = Vpairs[ii].first;
    int W = Vpairs[ii].second;
    
    pair<int,int> vpair = make_pair( V, W );
    it = lower_bound( Vpairs_begin, Vpairs_end, vpair );
    bool bad1 = ( it == Vpairs_end || ! ( *it == vpair ) );

    pair<int,int> rcpair = make_pair( toVRC[W], toVRC[V] );
    itRC = lower_bound( Vpairs_begin, Vpairs_end, rcpair );
    bool bad2 = ( itRC == Vpairs_end || ! ( *itRC == rcpair ) );

    if ( bad1 || bad2 ) {
      n_bad++;
      continue;
    }
    
    int e_id = eids[ distance( Vpairs_begin, it ) ];
    int rce_id = eids[ distance( Vpairs_begin, itRC ) ];
    toERC[e_id] = rce_id;
    toERC[rce_id] = e_id;
  }
  
  // Done.
  if ( n_bad < 1 ) {
    if ( VertRC ) *VertRC = toVRC;
    if ( EdgeRC ) *EdgeRC = toERC;
  }
  
  return ( n_bad == 0 );

}

/**
 * PrintBackbone
 */
void PrintBackbone( const backbone &bbg,
		    const String &out_file,
		    const vecbvec *bases )
{
  ofstream out( out_file.c_str( ) );

  const int n_edges = bbg.EdgeObjectCount( );
  const int n_vertices = bbg.N( );

  vec<int> to_left;
  vec<int> to_right;
  bbg.ToLeft( to_left );
  bbg.ToRight( to_right );

  out << "ADJACENCIES\n" << endl;
  if ( n_edges < 1 ) 
    out << "no edges (adjacencies) found\n";
  for (int ii=0; ii<n_edges; ii++)
    out << "V" << to_left[ii] << "  ==>  " << "V" << to_right[ii] << "\n";
  out << endl;

  out << "VERTICES\n" << endl;
  for (int ii=0; ii<n_vertices; ii++) {
    const vec<int> &vs = bbg.Vert( ii );
    out << "V" << ii << "  (" << vs.size( ) << " vertices)";
    if ( vs.isize( ) < 1 ) out << "\n";
    else out << ":  ";
    for (int jj=0; jj<vs.isize( ); jj++) {
      int uid = vs[jj] / 2;
      out << vs[jj]
	  << ( bases ? " [" + ToString( (*bases)[uid].size( ) ) + " bp]" : "" )
	  << ( jj == vs.isize( ) - 1 ? "\n" : " -> " );
    }
  }
  out << endl;
  
  out.close( );
}

/**
 * BackboneLinearize
 *
 * Compactify linear sequences from a given backbone. If
 *
 *   V_1 --> V_2 --> ... --> V_N
 *
 * with no other edges from/to any of the intermediate V_i's, then the
 * sequence above is replaced by the single vertex
 *
 *   V_1
 *
 * where V_1 is obtained by concatenating the old V_1, ..., V_N (the
 * other vertices are resized to 0).
 */
void BackboneLinearize( const backbone &in_backbone,
			backbone &out_backbone )
  
{
  // Core info.
  int n_vert = in_backbone.N( );
  int n_edges = in_backbone.EdgeObjectCount( );

  vec<int> to_left;
  vec<int> to_right;
  in_backbone.ToLeft( to_left );
  in_backbone.ToRight( to_right );
  
  // Count number of edges going from and to each vertex.
  vec<int> n_edges_in( n_vert, 0 );
  vec<int> n_edges_out( n_vert, 0 );
  for (int eid=0; eid<n_edges; eid++) {
    n_edges_in[ to_right[eid] ]++;
    n_edges_out[ to_left[eid] ]++;
  }
  
  // Descriptors for output graph.
  vec< vec<int> > from( n_vert );
  vec< vec<int> > to( n_vert );
  vec< vec<int> > fromEO( n_vert );
  vec< vec<int> > toEO( n_vert );
  vec< vec<int> > vert( n_vert );
  vec<CProx> edges;
  
  // Initialize sources and vertices (initially, source[ii] = ii).
  vec<int> source;
  source.reserve( n_vert );
  for (int ii=0; ii<n_vert; ii++) {
    vert[ii] = in_backbone.Vert( ii );   // may be empty
    if ( vert[ii].size( ) > 0 ) source.push_back( ii );
    else source.push_back( -1 );
  }

  // Fill vertices
  vec<bool> to_collapse( n_edges, false );
  for (int eid=0; eid<n_edges; eid++) {
    int v = to_left[eid];
    int w = to_right[eid];
    if ( n_edges_out[v] != 1 || n_edges_in[w] != 1 ) continue;
    
    // Skip bi-directional edges (v -> w, and w -> v).
    const vec<int> &local = in_backbone.From( w );
    if ( binary_search( local.begin( ), local.end( ), v ) ) continue;
    
    // Circular structure.
    if ( w == source[v] ) continue;
    
    // Collapse edge.
    int lstart = source[v];
    while ( lstart > -1 && source[lstart] != lstart  ) lstart = source[lstart];

    source[w] = lstart;
    copy( vert[w].begin(), vert[w].end(), back_inserter(vert[lstart]) );
    vert[w].clear( );
  }

  // Fill all other data structures. Note
  for (int eid=0; eid<n_edges; eid++) {
    int v = source[ to_left[eid] ];
    while ( v > -1 && source[v] != v ) v = source[v];

    int w = source[ to_right[eid] ];
    while ( w > -1 && source[w] != w ) w = source[w];

    if ( v < 0 || w < 0 ) continue;   // do not connect to/from empty vertices
    if ( v == w ) continue;   // this cuts circular structures
    if ( vert[v].size( ) < 1 ) continue;
    if ( vert[w].size( ) < 1 ) continue;
    
    from[v].push_back( w );
    fromEO[v].push_back( edges.size( ) );
    to[w].push_back( v );
    toEO[w].push_back( edges.size( ) );
    edges.push_back( in_backbone.EdgeObject( eid ) );
  }
  
  // Sort froms and tos.
  for (int ii=0; ii<n_vert; ii++) {
    sort( from[ii].begin( ), from[ii].end( ) );
    sort( to[ii].begin( ), to[ii].end( ) );
    sort( fromEO[ii].begin( ), fromEO[ii].end( ) );
    sort( toEO[ii].begin( ), toEO[ii].end( ) );
  }

  // Generate graph.
  out_backbone.Initialize( from, to, vert, edges, toEO, fromEO );
  ForceAssert( IsSymmetric( out_backbone ) );
  
}

/**
 * AreConsistent
 *
 * Decide if a triangle
 *
 *         e01       e12
 *     v0  --->  v1  ---> v2
 *
 *              e02
 *     v0  -------------> v2
 *
 * is consistent, ie if the gap implied by e02 is consistent with the
 * gap implied by ( e01 + sequence at v1 + e12 ).
 */
bool AreConsistent( const backbone &in_backbone,
		    const vec<int> &to_left,
		    const vec<int> &to_right,
		    const int e01,
		    const int e12,
		    const int e02 )
{
  // Check input makes sense.
  ForceAssertEq( to_right[e01], to_left[e12] );
  ForceAssertEq( to_left[e01], to_left[e02] );
  ForceAssertEq( to_right[e12], to_right[e02] );
  
  
  // SANTEMP - TODO - SANTEMP

  //   ...
  
  // It is consistent.
  return true;
}			 

/**
 * BackboneTriangularize
 *
 * Remove edge v1 -> v3 if there are edges v1 -> v2 -> v3 which form a
 * consistent triangle with v1 -> v3. It returns the number of edges
 * removed.
 */
int BackboneTriangularize( const backbone &in_backbone,
			   backbone &out_backbone,
			   ostream *log )
  
{
  // Core info.
  int n_vertices = in_backbone.N( );
  int n_edges = in_backbone.EdgeObjectCount( );

  if ( log ) *log << Date( ) << " [BL]: building indices" << endl;
  vec<int> to_left;
  vec<int> to_right;
  in_backbone.ToLeft( to_left );
  in_backbone.ToRight( to_right );

  // Edges tagged for removal.
  vec<bool> killed( n_edges, false );
  
  // Involution maps.
  vec<int> VertRC;
  vec<int> EdgeRC;
  ForceAssert( IsSymmetric( in_backbone, &VertRC, &EdgeRC ) );

  // Loop over all sets ( v0 -> v1 -> v2 ), ( v0 -> v2 ).
  int dotter = 1000;
  if ( log ) *log << Date( ) << " [BT]: testing "
		  << n_vertices << " vertices (.="
		  << dotter << " vertices)"
		  << endl;
  for (int v0=0; v0<n_vertices; v0++) {
    if ( log && v0 % dotter == 0 ) Dot( *log, v0 / dotter );
    const vec<int> &edges_01 = in_backbone.FromEdgeObj( v0 );

    vec< pair<int,int> > v1s;   // v1, id of edge v0 -> v1
    v1s.reserve( edges_01.size( ) );
    int n_edges_12 = 0;
    for (int eid=0; eid<edges_01.isize( ); eid++) {
      int v1 = to_right[ edges_01[eid] ];
      v1s.push_back( make_pair( v1, edges_01[eid] ) );
      n_edges_12 += in_backbone.FromEdgeObj( v1 ).size( );
    }
    sort( v1s.begin( ), v1s.end( ) );

    vec< pair<int,int> > v2s;   // v2, id of edge v1 -> v2
    v2s.reserve( n_edges_12 );
    for (int ii=0; ii<v1s.isize( ); ii++) {
      int v1 = v1s[ii].first;
      const vec<int> &edges = in_backbone.FromEdgeObj( v1 );
      for (int jj=0; jj<edges.isize( ); jj++)
	v2s.push_back( make_pair( to_right[ edges[jj] ], edges[jj] ) );
    }
    
    vec< pair<int,int> >::const_iterator it1;
    vec< pair<int,int> >::const_iterator it2;
    for (int ii=0; ii<v2s.isize( ); ii++) {
      int v2 = v2s[ii].first;
      it2 = lower_bound( v1s.begin( ), v1s.end( ), make_pair( v2, 0 ) );
      if ( it2 == v1s.end( ) ) continue;
      if ( it2->first != v2 ) continue;

      int v1 = to_left[ v2s[ii].second ];
      it1 = lower_bound( v1s.begin( ), v1s.end( ), make_pair( v1, 0 ) );
      ForceAssert( it1 != v1s.end( ) );
      ForceAssert( it1->first == v1 );

      int e12 = v2s[ii].second;    // edge v1 -> v2
      int e01 = it1->second;       // edge v0 -> v1
      int e02 = it2->second;       // edge v0 -> v2
      if ( ! AreConsistent( in_backbone, to_left, to_right, e01, e12, e02 ) )
	continue;

      // Do not erase a v -> w edge if any of the edges is bi-directional.
      int vi = to_right[e01];
      int v = to_left[e02];
      int w = to_right[e02];
      const vec<int> &fromws = in_backbone.From( w );
      const vec<int> &fromvis = in_backbone.From( vi );
      if ( binary_search( fromws.begin( ), fromws.end( ), v ) ) continue;
      if ( binary_search( fromws.begin( ), fromws.end( ), vi ) ) continue;
      if ( binary_search( fromvis.begin( ), fromvis.end( ), v ) ) continue;
      
      // Kill edges (both edge and its rc).
      killed[e02] = true;
      killed[ EdgeRC[e02] ] = true;
    }
  }
  if ( log ) *log << endl;

  // Rehab some killed edges (and their rcs). In this case:
  //    v ---> wi (i=1, 2, 3)
  //    w1 ---> w2 ---> w3
  //    w3 ---> w1
  // all edges v ---> wi would be tagged for deletion, which is not
  // what we want (we rehab all v ---> wi).
  for (int vid=0; vid<n_vertices; vid++) {
    const vec<int> &eidsf = in_backbone.FromEdgeObj( vid );
    bool all_killed = true;
    for (int ii=0; ii<eidsf.isize( ); ii++) {
      if ( ! killed[ eidsf[ii] ] ) {
	all_killed = false;
	break;
      }
    }
    if ( all_killed ) {
      for (int ii=0; ii<eidsf.isize( ); ii++) {
	killed[ eidsf[ii] ] = false;
	killed[ EdgeRC[ eidsf[ii] ] ] = false;
      }
    }

    const vec<int> &eidst = in_backbone.ToEdgeObj( vid );
    all_killed = true;
    for (int ii=0; ii<eidst.isize( ); ii++) {
      if ( ! killed[ eidst[ii] ] ) {
	all_killed = false;
	break;
      }
    }
    if ( all_killed ) {
      for (int ii=0; ii<eidst.isize( ); ii++) {
	killed[ eidst[ii] ] = false;
	killed[ EdgeRC[ eidst[ii] ] ] = false;
      }
    }
  }
  
  // Count killed events, and sanity check,.
  int n_killed = 0;
  for (int ii=0; ii<n_edges; ii++) {
    if ( killed[ii] ) {
      ForceAssert( killed[ EdgeRC[ii] ] );
      n_killed++;
    }
  }
  
  // Generate output graph.
  vec< vec<int> > from( n_vertices );
  vec< vec<int> > to( n_vertices );
  vec< vec<int> > fromEO( n_vertices );
  vec< vec<int> > toEO( n_vertices );
  vec< vec<int> > vertices;
  vec<CProx> edges;
  
  vertices.reserve( n_vertices );
  for (int ii=0; ii<n_vertices; ii++)
    vertices.push_back( in_backbone.Vert( ii ) );

  for (int eid=0; eid<n_edges; eid++) {
    if ( killed[eid] ) continue;
    int v = to_left[eid];
    int w = to_right[eid];
    from[v].push_back( w );
    fromEO[v].push_back( edges.size( ) );
    to[w].push_back( v );
    toEO[w].push_back( edges.size( ) );
    edges.push_back( in_backbone.EdgeObject( eid ) );
  }
  
  for (int ii=0; ii<n_vertices; ii++) {
    sort( from[ii].begin( ), from[ii].end( ) );
    sort( to[ii].begin( ), to[ii].end( ) );
    sort( fromEO[ii].begin( ), fromEO[ii].end( ) );
    sort( toEO[ii].begin( ), toEO[ii].end( ) );
  }

  out_backbone.Initialize( from, to, vertices, edges, toEO, fromEO );

  // Sanity check.
  if ( log ) *log << Date( ) << " [BT]: testing IsSymmetric" << endl;
  ForceAssert( IsSymmetric( out_backbone ) );

  // Done.
  if ( log ) *log << Date( ) << " [BT]: done ("
		  << n_killed << " edges tagged for removal)"
		  << endl;
  return n_killed;
}



/**
 * Template instantiations
 */
#include "graph/DigraphTemplate.h"
template backbone::digraphVE();
template backbone::digraphVE( const proxgraph&,
			      const vec< vec<int> >& );
template void backbone::Initialize( const vec< vec<int> >&,
				    const vec< vec<int> >&,
				    const vec< vec<int> >&,
				    const vec<CProx>&,
				    const vec< vec<int> >&,
				    const vec< vec<int> >& );
template int backbone::N() const;
template const vec<int>& backbone::Vert( int ) const;
template void backbone::readBinary( BinaryReader& );
template void backbone::writeBinary( BinaryWriter& ) const;
