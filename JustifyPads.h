// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
// 

// Move pads around.

#ifndef JUSTIFY_PADS_H
#define JUSTIFY_PADS_H

#include "Alignment.h"
#include "Basevector.h"



/*
 * JustifyPadsRight
 *
 * Define two alignments between the same reads to be equivalent if they
 * match up modulo pads shifting, provided the bases on the shifted pads do
 * not change. For example:
 *
 * alignment 1:
 *       ......CAGTATTCTAAAAATCG.....
 *           ......ATTCT*AAAATCGTAC.....
 *
 * alignment 2:
 *       ......CAGTATTCTAAAAATCG.....
 *           ......ATTCTAA*AATCGTAC.....
 *
 * JustifyPadsRight will standardize the given alignment by pushing on
 * the right (but always preserving equivalence) all the gaps in the align.
 * Notice that whole gaps will not be broken (i.e. either it is possible
 * to justify the gap as a whole, or no action will be taken), but gaps
 * may be merged together.
 */
void JustifyPadsRight( align &al,
		       const basevector &bases1,
		       const basevector &bases2 );



#endif
