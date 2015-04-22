///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__ASSISTED__FILTER_ALIGNS_FOR_ASSISTING__H
#define PATHS__ASSISTED__FILTER_ALIGNS_FOR_ASSISTING__H

#include "Vec.h"
#include "lookup/LookAlign.h"

/**
 * AlignsGap
 *
 * It returns the gap size between left and right (or NaN if the two
 * belong to different targets).
 */
int AlignsGap( const look_align &left, const look_align &right );

/**
 * ChainsGap
 *
 * It returns the gap size between left and right chains (or NaN if
 * the two belong to different targets).
 */
int ChainsGap( const vec<int> &left_chain,
	       const vec<int> &right_chain,
	       const vec<look_align> &aligns );

/**
 * FilterBigGaps
 *
 *                     u3   u4  u5
 *                     ---- --- ------
 *       ----------- -------------------     -------- -------
 *       u1          u2                      u6       u7
 *
 * If the gap between u2 and u6 is within the MIN_GAP to MAX_GAP
 * range, then remove aligns of u3, u4, and u5, thereby eliminating
 * the large negative gap between u2 and u3.
 *
 * aligns: both input and output
 */
void FilterBigGaps( const int MIN_GAP,
		    const int MAX_GAP,
		    vec<look_align> &aligns,
		    ostream *log = 0 );

#endif
