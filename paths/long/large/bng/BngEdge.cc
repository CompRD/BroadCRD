///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "paths/long/large/bng/BngEdge.h"

typedef bng_edge F1;
template digraphE<F1>::digraphE(const vec<vec<int> >&, const vec<vec<int> >&,
     const vec<F1>&, const vec<vec<int> >&, const vec<vec<int> >&, const Bool);
template const F1& digraphE<F1>::EdgeObject(int) const;
template F1& digraphE<F1>::EdgeObjectMutable(int);
template void digraphE<F1>::DeleteEdges(const vec<int>&);
template void digraphE<F1>::ComponentsE( vec< vec<int> >& ) const;
template vec<int> digraphE<F1>::RemoveDeadEdgeObjects( );
template void digraphE<F1>::RemoveEdgelessVertices( );
template int digraphE<F1>::EdgeObjectIndexByIndexTo(int, int) const;
template int digraphE<F1>::EdgeObjectIndexByIndexFrom(int, int) const;
template void digraphE<F1>::TransferEdges(int, int, unsigned char);
template void digraphE<F1>::DeleteEdgeTo(int, int);
template void digraphE<F1>::DeleteEdgeFrom(int, int);
template void digraphE<F1>::Reverse( );

template void digraphE<F1>::SplayVertex(int);
template int digraphE<F1>::AddEdge(int, int, bng_edge const&);

template digraphE<F1> digraphE<F1>::Subgraph( const vec<int>& v ) const;

template Bool digraphE<F1>::EdgePaths( const int v, const int w, vec< vec<int> >&,
     const int, const int, const int ) const;

template void DistancesToEndArr( digraphE<F1> const& G, vec<int> const& edgeLens,
                            int const max_dist, Bool const fw, vec<int>& D );

template void digraphE<bng_edge>::writeBinary(BinaryWriter&) const;
template void digraphE<bng_edge>::readBinary(BinaryReader&);

void bng_edge::writeBinary( BinaryWriter& writer ) const
{    writer.write(len_);
     writer.write(id_pos_);    }

void bng_edge::readBinary( BinaryReader& reader )
{    reader.read(&len_);
     reader.read(&id_pos_);    }
