///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////



#include "lookup/LookAlignCleanup.h"


int RemoveIfNotOverlapping (vec<look_align_plus> & aligns, 
			    vec<vec<int> > & alignIndices,
			    const vec<HoIntervalWithId> & valid){
  //check that the intervals are ordered and disjoint
  for (int i=0; i < valid.isize()-1; ++i) {
    Assert(LessById(valid[i], valid[i+1]) &&
           ( valid[i].id != valid[i+1].id || !Overlap(valid[i], valid[i+1])) );
  }
  //PRINT(valid);
  int removedCount=0;
  vec<HoIntervalWithId>::const_iterator vIter;
  for (size_t i=0; i != alignIndices.size(); ++i) {
    if (alignIndices[i].empty()) continue;
    const look_align & la = aligns[alignIndices[i][0]];
    HoIntervalWithId tofind(la.a.pos2(), la.a.Pos2(), la.target_id);
    vIter = upper_bound(valid.begin(), valid.end(), tofind, LessById);
    //PRINT4(tofind, *validIter, *valid.begin(), *valid.end());
    if (vIter == valid.begin() || !(--vIter)->Contains(tofind) ) {
      //remove all alignments for this read from the indices
      //because the best one is not in the valid intervals.
      alignIndices[i].clear();
      ++removedCount;
    }
  }
  return removedCount;
}


/// Use RemoveIfNotOverlapping according to the CALL_ONLY parameter.
int ProcessCallOnly(int i, const vec<String> & CALL_ONLY,
		    vec<look_align_plus> & aligns, 
		    vec<vec<int> > & alignIndices) {
  int removed = 0;
  if (CALL_ONLY.isize() > i  && IsRegularFile(CALL_ONLY[i])) {
    cout << "Reading CALL_ONLY file and cleaning alignments... " << endl;
    Ifstream(is, CALL_ONLY[i]); 
    vec<HoIntervalWithId> callOnly;
    callOnly.ReadFromTextStream(is);
    sort(callOnly.begin(), callOnly.end(), LessById);
    removed = RemoveIfNotOverlapping(aligns, alignIndices, callOnly);
    cout << removed << " reads removed because the best alignment "
	 << " was not in a high quality region of the reference." << endl;
  }
  return removed;
}

