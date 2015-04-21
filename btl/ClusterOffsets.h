///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef BTL__CLUSTER_OFFSETS__H
#define BTL__CLUSTER_OFFSETS__H

#include "Vec.h"
#include "math/Functions.h"

/**
 * ClusterOffsets
 *
 * Build clusters of offsets from a sorted vector of input offsets (as
 * generated, for example, from FindPossibleAlignments, in
 * pairwise_aligners/KmerAligner).
 *
 * MIN_SEEDS: discard clusters with less than MIN_SEEDS offsets in them
 * offsets: code will ForceAssert if this is not sorted
 * clusters: output (may be empty)
 */
void ClusterOffsets( const int MIN_SEEDS,
		     const vec<int> &offsets,
		     vec<int> &clusters );

#endif
