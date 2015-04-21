// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#include "pairwise_aligners/PatchAligner.h"



//==================================================================
void CAlignPatch::AsString(String & s) const
{
  unsigned int i;
  char tmp[16];
  tmp[1] = 0;
  for (i=0; i<m_seq.size(); i++) {
    tmp[0] = as_base(m_seq[i]);
    s += tmp;
  }
}

void CAlignPatch::Set(int start, int stop)
{
  m_start = start;
  m_stop = stop;
}


void CAlignPatch::SetSeq(int start, int stop, const basevector & b)
{
  //m_seq.clear();
  
  if (stop - start > 0)
    m_seq.SetToSubOf(b, start, stop-start);

  if (m_stop <= m_start) {
    start -= (m_start - m_stop) + 1;
    if (start < 0)
      start = 0;
    m_seq.SetToSubOf(b, start, stop-start);   
    m_start = m_stop - 1;
    //m_stop;
  }

}







void PrintPatchAlign(const vec<CAlignPatch> & patches, 
		     const vec<int> & matchPositions,
		     const basevector & seq,
		     const basevector & ref)
{

  int i, j, k;
  cout << "-------------- <alignment> ---------------------" << endl;
  for (i=0; i<matchPositions.isize(); i++) {
    if (matchPositions[i] != -1) {
      cout << i << " -> " << matchPositions[i];
      cout << ", " << as_base(seq[i]);
      cout << " -> " << as_base(ref[matchPositions[i]]); 
    }
    for (j=0; j<patches.isize(); j++) {
      if (patches[j].GetStart() == matchPositions[i]) {
	String alt;
	patches[j].AsString(alt);
	for (k=patches[j].GetStart(); k<patches[j].GetStop(); k++) {
	  cout << as_base(ref[k]);
	}
	cout << " ====> " << alt << endl;
      }
    }
    cout << endl;

  }
  cout << "-------------- </alignment> ---------------------" << endl << endl;

}
