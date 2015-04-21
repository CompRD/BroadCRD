///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef CLEANSE_SCAFFOLD_GRAPH_H
#define CLEANSE_SCAFFOLD_GRAPH_H

#include "Superb.h"
#include "Vec.h"
#include "graph/Digraph.h"
#include "paths/Alignlet.h"
#include "paths/BuildScaffoldLocs.h"
#include "paths/ContigsManager.h"
#include "paths/reporting/CLinkBundle.h"

/**
 * CleanseScaffoldGraph
 *
 * Iterative procedure to break bad supers. Look at all the supers
 * linked directly to a given super (the "center"), and decide based
 * on the spread of the links if and where to break the center. Data
 * structures will be updated. It returns the number of cutting
 * performed.
 *
 * LIGHT: if true only isolate supers (whitout setting index of aligns to -1)
 * MIN_LINKS: to connect a super to the center super
 * MAX_OVERLAP: tag a situation as conflict if there is a large implied overlap
 */
int CleanseScaffoldGraph( digraphE<CLinkBundle> &graph,
			  vec<superb> &supers,
			  ContigsManager &manager,
			  ostream &out,
			  const bool VERBOSE,
			  const bool LIGHT = False,
			  const int MIN_LINKS = 48,
			  const int MAX_OVERLAP = 350 );

#endif
