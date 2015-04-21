///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "lookup/LookAlignSort.h"
#include "lookup/LookAlign.h"
#include "Vec.h"
#include "lookup/FlowAlignSummary.h"

int GetReverseIndex(const vec<int> & queryIds, vec<int> & rindex, int offset)
{
  const int RINDEX_SIZE = *max_element(queryIds.begin(),queryIds.end()) + 1;
  rindex.assign(RINDEX_SIZE, -1);

  //figure out which index corresponds to which name
  for (int i=0; i != queryIds.isize(); ++i) {
    int id = queryIds[i];
    rindex[id] = i+offset;
  }
  return RINDEX_SIZE;
}

void GetSortedAlignIndices (const String & fname, 
			    vec<look_align_plus> & aligns, 
			    vec<vec<int> > & alignIndices,
			    vec<int> * readIds,
			    bool sortByMutation) {
  if (fname.Contains("qltout")) {
    LoadLookAlignPlus(fname, aligns);
  }
  else {
    vec<look_align> temp;
    LoadLookAlignBinary(fname, temp);
    aligns.assign(temp.begin(), temp.end());
  }

  if (0 == readIds) {
    int maxid = 0;
    for (int i=0; i != aligns.isize(); ++i) {
      maxid = max(maxid, aligns[i].query_id);
    }
    readIds = new vec<int>(maxid+1,vec<int>::IDENTITY);
  }

  GetAlignIndices(*readIds, aligns, alignIndices);
  if (sortByMutation) {
    SortAlignIndices(aligns, alignIndices);
  }
  else {
    SortAlignIndicesByErrors(aligns, alignIndices);
  }
}

// TODO: potentially dangerous truncation of index by read
bool ApproveUniqueLookAligns( const vec<vec<int> > & alignIndices, 
			      const vec<look_align_plus> & aligns, int read,
			      double ALIGN_ERR_RATE, double ALIGN_COMP,
			      int * unalignedReads, 
			      int * badReads, 
			      int * ambigReads) { 
     
  if (alignIndices[read].empty()) {
    if (unalignedReads) ++*unalignedReads;
    return false;
  }
  const look_align & la = aligns[alignIndices[read][0]];
  if (la.ErrorRate() > ALIGN_ERR_RATE) {
    if (badReads) ++*badReads;
    return false;
  }
  if (alignIndices[read].size() > 1 
      && IsAmbiguous(la, aligns[alignIndices[read][1]], 0, ALIGN_COMP, True) ) {
    if (ambigReads) ++*ambigReads;
    return false; 
  }
  return true;
}

void GetBestAlignIndices (const String & fname, 
			  vec<look_align_plus> & aligns, 
			  vec< int > & bestIndices,
			  vec<int> * readIds,
			  bool sortByMutation) {
  if (fname.Contains("qltout")) {
    LoadLookAlignPlus(fname, aligns);
  }
  else {
    vec<look_align> temp;
    LoadLookAlignBinary(fname, temp);
    aligns.assign(temp.begin(), temp.end());
  }

  if (0 == readIds) {
    int maxid = 0;
    for (int i=0; i != aligns.isize(); ++i) {
      maxid = max(maxid, aligns[i].query_id);
    }
    readIds = new vec<int>(maxid+1,vec<int>::IDENTITY);
  }

  vec<vec<int> > alignIndices;
  GetAlignIndices(*readIds, aligns, alignIndices);
  if (sortByMutation) {
    BestAlignIndices(aligns, alignIndices, bestIndices);
  }
  else {
    BestAlignIndicesByErrors(aligns, alignIndices, bestIndices);
  }
}




