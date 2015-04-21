/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef MERGE_SCAFFOLD_GRAPH_H
#define MERGE_SCAFFOLD_GRAPH_H

#include "PairsManager.h"
#include "Superb.h"
#include "paths/Alignlet.h"
#include "paths/reporting/CLinkBundle.h"

/**
 * MergeScaffoldGraph
 *
 * Align edges which appear to overlap (based on linking evidence),
 * and if they overlap merge them. Both contigs and aligns are
 * changed.
 *
 * MIN_OVERLAP: only align edges if linking-based overlap is >= MIN_OVERLAP
 */
void MergeScaffoldGraph( vec<superb> &supers,
			 vec<fastavector> &contigs,
			 vec<alignlet> &aligns,
			 vec<int> &index,
			 ostream &out,
			 const PairsManager &pairs,
			 const bool VERBOSE = False,
			 const int MIN_OVERLAP = 500 );

#endif
