///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LOOK_ALIGN_CLEANUP_H
#define LOOK_ALIGN_CLEANUP_H

/** These methods clean up a vec<look_align_plus> and vec<vec<int> > of indices.

\file LookAlignCleanup.h

The vecs are expected to come from GetSortedAlignIndices or some related method
in LookAlignSort.h

*/

#include "lookup/LookAlign.h"
#include "Vec.h"

#include <functional>

/// Remove all indices if the best alignment is not in the valid intervals.
/// Return the number of reads for which indices were cleared.
int RemoveIfNotOverlapping (vec<look_align_plus> & aligns, 
			    vec<vec<int> > & alignIndices,
			    const vec<HoIntervalWithId> & valid);

/// Use RemoveIfNotOverlapping according to the CALL_ONLY parameter at index i.
/// The CALL_ONLY list may be incomplete or have fake file names,
/// such as "none", in which case we just do nothing at that index.
int ProcessCallOnly(int i, const vec<String> & CALL_ONLY,
		    vec<look_align_plus> & aligns, 
		    vec<vec<int> > & alignIndices);

/// Return true if la does not overlap valid.
/// Assumes all aligns come from the same read. 
inline bool RejectAlign(const vec<HoIntervalWithId> & valid, look_align & la) {
  HoIntervalWithId tofind(la.a.pos2(), la.a.Pos2(), la.target_id);
  vec<HoIntervalWithId>::const_iterator vIter = 
    upper_bound(valid.begin(), valid.end(), tofind, LessById);
    //PRINT4(tofind, *validIter, *valid.begin(), *valid.end());
  return (vIter == valid.begin() || !(--vIter)->Contains(tofind) );
}


#endif // LOOK_ALIGN_CLEANUP_H
