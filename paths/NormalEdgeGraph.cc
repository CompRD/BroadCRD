/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "paths/NormalEdgeGraph.h"

//
// Method sof class NormalEdgeGraph
//

NormalEdgeGraph::NormalEdgeGraph( const HyperKmerPath& hkp, nkmers_t MIN_NORMAL_EDGE_SIZE ) {

  minNormalEdgeSize = MIN_NORMAL_EDGE_SIZE;
  
  const vec< KmerPath >& edges = hkp.Edges();
  edge2normal.resize( edges.size(), -1 );
  for ( int e = 0; e < edges.isize(); e++ )
    if ( edges[ e ].KmerCount() > MIN_NORMAL_EDGE_SIZE ) {
      edge2normal[ e ] = normalEdges.isize();
      normalEdges.push_back( e );
    }

  vec<int> to_left, to_right;
  hkp.ToLeft( to_left );
  hkp.ToRight( to_right );

  vec< Bool > seenVerts( hkp.N() );
  matrix< Bool > normalEdgeAdjs( normalEdges.size(), normalEdges.size(), False );
  PRINT( seenVerts.size() );
  for ( normal_edge_id_t ne = 0; ne < normalEdges.isize(); ne++ ) {

    PRINT2( ne, normalEdges[ ne ] );

    fill( seenVerts.begin(), seenVerts.end(), False );
    int vertexAtTipOfNe = to_right[ normalEdges[ ne ] ] ;
    vec< int > activeVerts( 1, vertexAtTipOfNe );
    seenVerts[ vertexAtTipOfNe ] = True;

    while( activeVerts.nonempty() ) {
      int v = activeVerts.back();
      PRINT2( activeVerts.size(), v );
      activeVerts.pop_back();
      for ( int i = 0; i < hkp.From( v ).isize(); i++ ) {
	int w = hkp.From( v )[i];
	int e = hkp.EdgeObjectIndexByIndexFrom( v, i );
	ForceAssertEq( to_left[ e ], v );
	ForceAssertEq( to_right[ e ], w );
	if ( edge2normal[ e ] != -1 ) {
	  normalEdgeAdjs( ne, edge2normal[ e ] ) = True;
	} else {
	  if ( !seenVerts[ w ] ) {
	    activeVerts.push_back( w );
	    seenVerts[ w ] = True;
	  }
	}
      }
    }
    
  }

  Initialize( normalEdgeAdjs );
}

normal_edge_id_t NormalEdgeGraph::NormalEdgeId( int e ) const {
  return edge2normal[ e ];
}

void NormalEdgeGraph::writeBinary( BinaryWriter& writer ) const
{
    writer.write(minNormalEdgeSize);
    writer.write(normalEdges);
    writer.write(edge2normal);
    writer.write(*static_cast<digraph const*>(this));
}

void NormalEdgeGraph::readBinary( BinaryReader& reader )
{
    reader.read(&minNormalEdgeSize);
    reader.read(&normalEdges);
    reader.read(&edge2normal);
    reader.read(static_cast<digraph*>(this));
}
