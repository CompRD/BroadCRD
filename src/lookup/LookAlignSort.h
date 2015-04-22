///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LOOK_ALIGN_SORT_H
#define LOOK_ALIGN_SORT_H

#include "lookup/LookAlign.h"
#include "Vec.h"
#include "VecString.h"

#include <functional>

/// Read in look_aligns from file, sort, and load into vecs
/// \param fname name of look_align text file
/// \param aligns vector into which to put the alignments
/// \param alignIndices vec of sorted indices produced by GetAlignIndices
/// and SortAlignIndices.
/// \param readIds vector containing the ids of the reads in the order in
/// which they should go into alignIndices. That is, if readIds contains
/// {5, 7, 9}, the alignIndices will have size 3, and the first vec will
/// contain the indices in aligns of all the alignments for which query_id
/// is equal to 5. If the pointer is null, then it is assumed that the
/// vector is intended to be {0, 1, 2, ... n} where n is the largest
/// query_id in the alignments.

void GetSortedAlignIndices (const String & fname, 
			    vec<look_align_plus> & aligns, 
			    vec<vec<int> > & alignIndices,
			    vec<int> * readIds = 0,
			    bool sortByMutation = true);

/// Return true if there is a unique approved alignment for this read.
/// Preconditions: aligns and alignIndices have been created using
/// GetSortedAlignIndices with sortByMutation = false.
/// If the xxxReads pointers are not null, their pointees will be
/// incremented when we find reads in those categories.
/// "unique" and "approved are defined based on ALIGN_COMP and ALIGN_ERR_RATE.
bool ApproveUniqueLookAligns( const vec<vec<int> > & alignIndices, 
			      const vec<look_align_plus> & aligns, int read,
			      double ALIGN_ERR_RATE, double ALIGN_COMP,
			      int * unalignedReads = 0, 
			      int * badReads = 0, 
			      int * ambigReads = 0);


/// Read in look_aligns from file, find best, and load into vecs
/// \param fname name of look_align text file
/// \param aligns vector into which to put the alignments
/// \param alignIndices vec of indices to the best alignment for each read.
/// \param readIds vector containing the ids of the reads in the order in
/// which they should go into alignIndices. That is, if readIds contains
/// {5, 7, 9}, the alignIndices will have size 3, and the first vec will
/// contain the indices in aligns of all the alignments for which query_id
/// is equal to 5. If the pointer is null, then it is assumed that the
/// vector is intended to be {0, 1, 2, ... n} where n is the largest
/// query_id in the alignments.

void GetBestAlignIndices (const String & fname, 
			    vec<look_align_plus> & aligns, 
			    vec< int > & bestIndices,
			    vec<int> * readIds = 0,
			    bool sortByMutation = true);



/// Given a list of queryIDs, queryIDs, compute the reverse index,
/// rindex, from the queryId to the position in the list.  The reverse
/// index is defined by the property rindex[queryIDs[i]]==i if
/// 0<=i<queryIDs.size().  In addition queryIDs[rindex[i]]==i if
/// rindex[i]>=0.  The return value is the largest queryID found,
/// which becomes the size of rindex.  The optional offset is added to
/// all values in rindex.
/// If a value in rindex would not normally be filled, it is set to -1.
int GetReverseIndex(const vec<int> & queryIds, vec<int> & rindex,int offset=0);
///Sort and uniqify all query ids in the align vec into a vec.
template <class AlignVec>
void GetSortedQueryIds(const AlignVec & aligns, vec<int> & v) {
  set<int> ids;
  for (int i=0; i != aligns.isize(); ++i) {
    ids.insert(aligns[i].query_id);
  }
  v.assign(ids.begin(), ids.end());
}

/// Get all the alignments and put their indices into a vector organized 
/// by the query_id, as contained in the look_align and in names.
/// Note that alignIndices is indexed into by id, that is by position in
/// the vector queryIds, and not by the int contained at that position.
template <class AlignVec>
void GetAlignIndices(const vec<int> & queryIds, 
		     const AlignVec & aligns,
		     vec<vec<int> > & alignIndices) {
  vec<int> rindex;
  const int RINDEX_SIZE = GetReverseIndex(queryIds, rindex);
  const int N = queryIds.size();

  //save the aligns for each name in the vec for each index.
  alignIndices.assign(N,vec<int>());
  const int A=aligns.size();
  for (int a = 0; a !=A; ++a) {
    int id=aligns[a].query_id;
    if (id >=0 && id < RINDEX_SIZE && rindex[id] != -1 ) { 
      //make sure this align is for something
      //that was in the queryIds vector!
      alignIndices[rindex[id]].push_back(a);
    }
  }
}

/// Get all the alignments and put their indices into a vector organized 
/// by the query_id, as contained in the look_align and in names.
/// Note that alignIndices is indexed into by id, that is by position in
/// the vector names, and not by the int contained at that position in names.
template <class AlignVec>
void GetAlignIndices(vecString & names, 
		      const AlignVec & aligns,
		      vec<vec<int> > & alignIndices) {
  vec<int> inames;
  transform(names.begin(), names.end(), back_inserter(inames),
	    mem_fun_ref(&String::Int));
  GetAlignIndices(inames, aligns, alignIndices);
}

/// Struct for sorting a vec of look_align or look_align_plus for the same query.
/// Sort is based on mutation rate, but stable for equal mutation rates.
template <class AlignVec>
struct BetterLookAlign: public std::binary_function<int,int,bool> {
  BetterLookAlign(const AlignVec & aligns): aligns(aligns) {}
  bool operator()(int l, int r) {
    double le = aligns[l].MutationRate();
    double re = aligns[r].MutationRate();
    if (le < re) return true;
    if (re < le) return false;
    return (l<r);
  }
private:
  const AlignVec & aligns;
};

/// Sort by error rate, not mutation rate.
template <class AlignVec>
struct SortLookAlignByErrors: public std::binary_function<int,int,bool> {
  SortLookAlignByErrors(const AlignVec & aligns): aligns(aligns) {}
  bool operator()(int l, int r) {
    double le = aligns[l].ErrorRate();
    double re = aligns[r].ErrorRate();
    if (le < re) return true;
    if (re < le) return false;
    return (l<r);
  }
private:
  const AlignVec & aligns;
};

template <class AlignVec, class Comp>
void SortAlignIndicesInternal(const AlignVec & aligns, 
			      vec<vec<int> > & alignIndices, Comp comp) {
  const int N = alignIndices.isize();
  for (int i=0; i != N; ++i) {
    sort(alignIndices[i].begin(), alignIndices[i].end(), comp);
  }
}  

/// Sort vectors so the best alignment for each read is at alignIndices[i][0].
/// Best is defined as having the lowest MutationRate().
template <class AlignVec >
void SortAlignIndices(const AlignVec & aligns, vec<vec<int> > & alignIndices) {
  BetterLookAlign<AlignVec> comp(aligns);
  SortAlignIndicesInternal(aligns, alignIndices, comp);
}

template <class AlignVec >
void SortAlignIndicesByErrors(const AlignVec & aligns, vec<vec<int> > & alignIndices) {
  SortLookAlignByErrors<AlignVec> comp(aligns);
  SortAlignIndicesInternal(aligns, alignIndices, comp);
}


template <class AlignVec, class Comp >
void BestAlignIndicesInternal(const AlignVec & aligns, 
			      vec<vec<int> > & alignIndices, 
			      vec<int> & bestIndices, Comp comp) {
  const int N = alignIndices.isize();
  bestIndices.assign(N,-1);
  for (int i=0; i != N; ++i) {
    if ( alignIndices[i].empty() ) continue;
    bestIndices[i] = *min_element(alignIndices[i].begin(), 
				  alignIndices[i].end(), comp);
  }
}

/// Sort vectors so the best alignment for each read is at alignIndices[i][0].
/// Best is defined as having the highest score() when turned into a 
/// FlowAlignSummary.
template <class AlignVec >
void BestAlignIndices(const AlignVec & aligns, vec<vec<int> > & alignIndices, vec<int> & bestIndices) {
  BetterLookAlign<AlignVec> comp(aligns);
  BestAlignIndicesInternal(aligns, alignIndices, bestIndices, comp);
}


template <class AlignVec >
void BestAlignIndicesByErrors(const AlignVec & aligns, vec<vec<int> > & alignIndices, vec<int> & bestIndices) {
  SortLookAlignByErrors<AlignVec> comp(aligns);
  BestAlignIndicesInternal(aligns, alignIndices, bestIndices, comp);
}


#endif //LOOK_ALIGN_SORT_H
