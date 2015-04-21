/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef REMOVE_LOCS_OVER_WINDOW_H
#define REMOVE_LOCS_OVER_WINDOW_H

#include "LocsHandler.h"
#include "SeqInterval.h"
#include "Vec.h"

/**
 * RemoveLocsOverWindows
 *
 * Remove locs from a locs file. Argument algo_select is used to decide
 * the criterion for removal.
 *
 * algo_select = 1: remove all locs that overlap any of the given windows. 
 *
 * algo_select: see above
 * wins: windows on contigs (seq_id = contig_id)
 * locs_in: input locs
 * locs_out: output locs
 * plog: optional log stream
 */
void RemoveLocsOverWindows( const int algo_select,
			    const vec<seq_interval> &wins,
			    const lhandler &locs_in,
			    vec<read_location> &locs_out,
			    ostream *plog = 0 );

#endif
