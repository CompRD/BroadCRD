// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#ifndef _PATCHALIGNER_H_
#define _PATCHALIGNER_H_


#include "Vec.h"
#include "String.h"
#include "Basevector.h"
#include "pairwise_aligners/KmerAligner.h"

//#define PATCH_ALIGNER_K 48

//------------------------------------------------------------
class CAlignPatch
{
 public:
  CAlignPatch() {
    m_start = m_stop = -1;
  }

  ~CAlignPatch() {}


  int GetStart() const {return m_start;}
  int GetStop() const {return m_stop;}

  const basevector & GetSequence() const {return m_seq;}
  void AsString(String & s) const;


  void Set(int start, int stop);
  // This is start & stop on the SEQUENCE!
  void SetSeq(int start, int stop, const basevector & b);

  bool operator == (const CAlignPatch & patch) const {
    if (m_start == patch.m_start 
	&& m_stop == patch.m_stop
	&& m_seq == patch.m_seq)
      return true;
    else
      return false;
  }

  bool operator < (const CAlignPatch & patch) const {
    if (m_start < patch.m_start)
      return true;
    if (m_start > patch.m_start)
      return false;
    if (m_stop < patch.m_stop)
      return true;
    if (m_stop > patch.m_stop)
      return false;
    return (m_seq < patch.m_seq);
  }

 private:
  int m_start;
  int m_stop;
  basevector m_seq;
  
};




//------------------------------------------------------------
template <int PATCH_ALIGNER_K=48>
class CPatchAligner
{
 public:
  CPatchAligner() {
    m_slack = 48;
  }

  void SetSlack(int i) {m_slack = i;}

  ~CPatchAligner() {}

   void SetBases(const basevector & ref) {
     m_kmerAligner.SetBases(ref);
   }


   bool Align(vec<CAlignPatch> & patches, 
	      vec<int> & matchPositions,
	      const basevector & seq,
	      const basevector & ref,
	      int suspectedPos);

 private:
  int GetClosest(const basevector & kmerbases,
		 int pos);



  int m_slack;
  KmerAligner<PATCH_ALIGNER_K> m_kmerAligner;

};



//==================================================================
template <int PATCH_ALIGNER_K>
bool CPatchAligner<PATCH_ALIGNER_K>::Align(vec<CAlignPatch> & patches, 
					   vec<int> & matchPositions,
					   const basevector & seq,
					   const basevector & ref,
					   int suspectedPos)
{
  //KmerAligner<PATCH_ALIGNER_K> kmerAligner;

  //kmerAligner.SetBases(ref);

  vec<int> positions;
  vec<int> ids;
	

  
  int i, j;

  basevector kmerbases;
  matchPositions.clear();
  matchPositions.resize(seq.size(), -1);


  int patchStart = 0;
  int patchEnd = (int)seq.size();
  int patchStartRef = suspectedPos;
  int patchEndRef = suspectedPos + (int)seq.size();

  bool bNewStart = true;

  for (i=0; i<(int)seq.size() - PATCH_ALIGNER_K; i++) {
    kmerbases.SetToSubOf(seq, i, PATCH_ALIGNER_K);

    int pos = GetClosest(kmerbases,
			 suspectedPos + i);
    // Found a match...
    if (pos != -1) {
      patchEnd = i;
      patchEndRef = pos;

      // Adjust position for hanging ends...
      if (patchStart == 0) {
	patchStartRef = pos - patchEnd;
      }
      
      if (patchEndRef != patchStartRef || patchStart != patchEnd) {
	CAlignPatch patch;
	patch.Set(patchStartRef, patchEndRef);
	patch.SetSeq(patchStart, patchEnd, seq);
	patches.push_back(patch);
      }

      bNewStart = true;

      //for (j=i; j<i+PATCH_ALIGNER_K; j++)
      //matchPositions[j] = pos + j;

      //i += PATCH_ALIGNER_K;
      //pos += PATCH_ALIGNER_K;
      while (i<(int)seq.size() && pos <(int)ref.size()) {
	if (seq[i] != ref[pos])
	  break;

	//cout << "assigning " << i << " with  " << pos << endl;
	matchPositions[i] = pos;
	i++;
	pos++;
      }
    }
    if (bNewStart) {
      patchStart = i;
      patchStartRef = pos;
      bNewStart = false;
    }
  }

  // Hanging end??
  if (patchStart + 1 < (int)seq.size()) {
    patchEnd = (int)seq.size();
    patchEndRef = patchStartRef + patchEnd - patchStart;
    CAlignPatch patch;
    patch.Set(patchStartRef, patchEndRef);
    patch.SetSeq(patchStart, patchEnd, seq);
    patches.push_back(patch);
  }

  return true;

}


template <int PATCH_ALIGNER_K>
int CPatchAligner<PATCH_ALIGNER_K>::GetClosest(const basevector & kmerbases,
					       int pos)
{


  vec<int> ids;
  vec<int> positions;

  m_kmerAligner.FindPossibleAlignments(kmerbases, 
				       positions, 
				       ids);

  int i;
  int closestDiff = pos + 1000 * m_slack;
  int closestPos = -1;
 
  for (i=0; i<positions.isize(); i++) {
    int diff = positions[i] - pos;
    if (diff < 0)
      diff = -diff;
    if (diff > m_slack)
      continue;
    if (diff < closestDiff) {
      closestDiff = diff;
      closestPos = positions[i];
    }
  }
  return closestPos;
}




void PrintPatchAlign(const vec<CAlignPatch> & patches, 
		     const vec<int> & matchPositions,
		     const basevector & seq,
		     const basevector & ref);


#endif //_PATCHALIGNER_H_

