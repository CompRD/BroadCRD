// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
//

#ifndef DETECT_VALID_PAIRS_H
#define DETECT_VALID_PAIRS_H

#include "ReadPairing.h"
#include "lookup/LookAlign.h"

/**
 * DetectValidPairs
 *
 * Find all possible valid pairs. Notice that reads may own multiple aligns,
 * in which case more than one valid pairs may be returned. Warning: hits
 * must be sorted by query_id (or it will ForceAssert).
 *
 * If both max_stretch and max_mult are < 0, then any stretch will be
 * allowed. This amounts to ask for all the logical pairs.
 *
 * to_pair: map sending a read_id to the id of the pair containing it
 * valid: output (positions in hits of the two end reads of a valid pair)
 * max_stretch: max stretch (argument of IsValidPair( ))
 * max_mult: max mult stretch (argument of IsValidPair( ) ))
 * cap: cap for valid pairs to be found for a given pair id
 */
void DetectValidPairs( const vec<look_align_plus> &hits,
		       const vec<read_pairing> &pairs,
		       const vec<int> &to_pair,
		       vec< pair<int,int> > &valid,
		       double max_stretch = 4.5,
		       double max_mult = 2.0,
		       int *cap = 0 );

#endif
