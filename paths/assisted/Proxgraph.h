///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__ASSISTED__PROXGRAPH__H
#define PATHS__ASSISTED__PROXGRAPH__H

#include "graph/Digraph.h"
#include "paths/assisted/CProx.h"

/**
 * proxgraph
 *
 * A proximity graph. Please, only template instantiations and high
 * level functions, in here!
 */
typedef digraphE<CProx> proxgraph;
extern template class digraphE<CProx>;

/**
 * IsSymmetric
 *
 * Check that a proxgraph is symmetric, ie an edge e: v -> w exists if
 * and only if the edge e': rc(w) -> rc(v) exists.
 */
bool IsSymmetric( const proxgraph &pgraph );

/**
 * AreIdentical
 *
 * Check if two proxgraphs are identical (if log is not null, send
 * progress to log).
 */
bool AreIdentical( const proxgraph &left,
		   const proxgraph &right,
		   ostream *log = 0 );

/**
 * RemoveContraryRefLinks: remove edges from verticies which are
 * supported by reference only where there are other edges supported
 * by links.
 */

void RemoveContraryRefLinks(proxgraph &pgraph);
#endif
