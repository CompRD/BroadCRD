///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef BTL__WALK_ALIGN__H
#define BTL__WALK_ALIGN__H

#include "Alignment.h"

/**
 * WalkOn1
 *
 * Given an align al between read1 and read2, find pos on2 on read2
 * corresponding to given pos on1 on read1. Note that class align
 * already has methods PosOn1 and PosOn2, but these do not deal with
 * gaps.
 *
 * WalkOn1 returns a pair of ints, where "first" captures on2, and
 * "second" the distance in the gap between where on1 lies on the gap
 * and on2 (hence, it is 0 unless pos1 corresponds to a gap in
 * read2). For example:
 *
 *        0123456789...              on1
 *        ACGGTACG TT ACTATTT        read1
 *     AAAACGCT  GTTTTAGTAT          read2
 *     01234567  89...               on2
 *
 * align:  pos1  pos2  (gap,len)  (gap,len)  (gap,len)  (gap,len)
 *         0     3     (0,5)      (-2,1)     (1,2)      (1,5)
 *
 * WalkOn1 returns:
 *
 *    on1  ->  <on2,extra>
 *      0      <3,0>
 *      4      <7,0>
 *      5      <8,2>
 *      6      <8,1>
 *      7      <8,0>
 *      8      <9,0>
 *
 * WalkOn1 returns -1 if on1 < pos1, and nan if on1 > Pos1.
 */
pair<int,int> WalkOn1( const align &al, const int on1 );


/**
 * WalkOn2
 *
 * Same as WalkOn1 (mutatis mutandis: replace on1 with on2)
 */
pair<int,int> WalkOn2( const align &al, const int on1 );

#endif
