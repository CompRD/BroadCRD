/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef VEC_TILING_H
#define VEC_TILING_H

#include "Vec.h"
#include "tiled/Tiling.h"

/**
 * LoadVecTiling
 *
 * Load a vector of tilings from file. If contig_ids is given, only load
 * tilings for the specified contigs (contig_ids must be sorted).
 */
void LoadVecTiling( const String &tilings_file,
		    vec<tiling> &tiles,
		    const vec<int> *contig_ids = 0 );

/**
 * SaveVecTiling
 *
 * Save a vector of tilings from file.
 */
void SaveVecTiling( const String &tilings_file,
		    const vec<tiling> &tiles );

#endif
