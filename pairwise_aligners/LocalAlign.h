// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#ifndef LOCAL_ALIGN_H
#define LOCAL_ALIGN_H

#include "Basevector.h"
#include "CoreTools.h"
#include "PackAlign.h"

/// LocalAlign( S, T, a, match, mismatch, gap )
///
/// Find the best local alignment of basevector "S" and "T" and put it in 
/// align "a".
///
/// This is computed relative to the given scores for "match", "mismatch",
/// and "gap", which default to 1, -1 and -2
///
/// If the best alignment has length 0, the alignment will have zero blocks, 
/// but be otherwise meaningless.
int LocalAlign( const basevector& S, const basevector& T, align& a,
		const int match = 1, 
		const int mismatch = -1, 
		const int gap = -2);

#endif
