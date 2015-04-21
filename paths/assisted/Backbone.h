///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__ASSISTED__BACKBONE__H
#define PATHS__ASSISTED__BACKBONE__H

#include "Basevector.h"
#include "graph/Digraph.h"
#include "paths/assisted/CProx.h"
#include "paths/assisted/Proxgraph.h"

/**
 * Backbone
 *
 * A backbone is a digraphVE in which vertices are lists of vertices
 * from some constant proxgraph (these can be thought of as
 * collapsed walks), and edges are CProx objects.
 *
 * A trivial backbone has the same number of vertices as the constant
 * proxgraph it is derived from, and v[ii][0] = { ii }, for all
 * ii's.
 *
 * Please: only template instantiations and high level functions, in
 * here!
 */
typedef digraphVE< vec<int>, CProx > backbone;
extern template class digraphVE<vec<int>,CProx>;

/**
 * BackboneFromCProx
 *
 * It takes in input a proxgraph and it generates in output a
 * trivial backbone digraphVE (see Backbone.h for details).
 */
void BackboneFromCProx( const proxgraph &in_graph,
			backbone &out_graph );

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
		      proxgraph &out_graphE );
/**
 * IsSymmetric
 *
 * Check that a backbone is symmetric:
 *
 *  1.  for each vertex V = {v1, ..., vn} there exists a unique vertex
 *      W != V such that W = rc(V), ie W = { rc(vn), ..., rc(v1) };
 *  2.  for each edge e: V -> W there exists a unique edge rc(e):
 *      rc(W) -> rc(V).
 *
 * If not null, fill VertRC (and/or EdgeRC) with the map V -> rc(V)
 * (with the map e -> rc(e)).
 */
bool IsSymmetric( const backbone &bbg,
		  vec<int> *VertRC = 0,
		  vec<int> *EdgeRC = 0 );

/**
 * PrintBackbone
 *
 * Human readable, text based.
 */
void PrintBackbone( const backbone &bbg,
		    const String &out_file,
		    const vecbvec *bases = 0 );

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
			backbone &out_backbone );


/**
 * BackboneTriangularize
 *
 * Remove edge v1 -> v3 if there are edges v1 -> v2 -> v3 which form a
 * consistent triangle with v1 -> v3. It returns the number of edges
 * removed.
 */
int BackboneTriangularize( const backbone &in_backbone,
			   backbone &out_backbone,
			   ostream *log = 0 );
#endif
