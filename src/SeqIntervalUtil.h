/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef SEQ_INTERVAL_UTIL_H
#define SEQ_INTERVAL_UTIL_H

#include "LocsHandler.h"
#include "SeqInterval.h"
#include "Vec.h"

/**
 * IntersectIntervalExternal
 *
 * Assume the given interval is an interval on a contig (so SeqId( )
 * is the contig_id). IntersectIntervalExternal returns the set of all
 * locs which have nonempty intersection with the interval.
 */
vec<int> IntersectIntervalExternal( const seq_interval &interval,
				    const lhandler &locs );

/**
 * IntersectIntervalInternal
 *
 * Same as IntersectIntervalExternal, the difference being that locs are
 * returned iff they are completely internal to the interval.
 */
vec<int> IntersectIntervalInternal( const seq_interval &interval,
				    const lhandler &locs );

/**
 * FlattenSeqIntervals
 *
 * Given a sorted VEC<seq_interval>, make a new one with overlapping
 * intervals merged -- so the same regions are covered, but all at
 * depth of coverage 1.  Each merged seq_interval gets the IntervalId
 * of the first seq_interval that contributed to it.
 */
void FlattenSeqIntervals( const vec<seq_interval>& to_flatten,
			  vec<seq_interval>& ans );

/**
 * SelectIntersecting
 *
 * Select from in_wins all and only the seq_intervals that intersect
 * at least one of the intervals in annotations. Save the output in
 * out_wins. Remark: annotations must be sorted.
 */
void SelectIntersecting( unsigned long n_contigs,
			 const vec<seq_interval>& in_wins,
			 const vec<seq_interval>& annotations,
			 vec<seq_interval> &out_wins );

/**
 * PrintSeqIntBasicNStats
 *
 * Run PrintBasicNStats on the lengths specified by the intervals.
 */
void PrintSeqIntBasicNStats( const vec<seq_interval> &intervals,
			     const String &name,
			     ostream &out );
  
#endif
