/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "pairwise_aligners/CombAligner.h"

int Combmer::m_k = 10;
int Combmer::m_skip = 2;
bool Combmer::m_bDouble = true;

void AddCombmers(vec<Combmer> & v, 
		 const basevector & ref,
		 int contigID)
{
  int k = Combmer::GetSize();

  int i, j;
  
  int start = v.isize();
  int newSize = start + (int)ref.size() - k + 2;
  if (newSize <= 0) {
    cout << "Error: trying to resize to " << newSize << ", ref size=" << ref.size() << endl;
    return;
  }

  v.resize(newSize);
  //cout << "ref size=" << ref.isize() << endl;
  //cout << "k=" << k << endl;

  for (i=0; i<(int)ref.size()-k+2; i++) {
    v[i+start].Set(ref, i, contigID, i);
  }  
}


int IsGoodEnough(const basevector & ref, 
		 const basevector & seq, 
		 int startOnRef,
		 bool bPrint, 
		 int & snps)
{
  //Do not allow for partial overlaps for now.
  if (startOnRef < 1)
    return false;
  if (startOnRef + (int)seq.size() > (int)ref.size())
    return false;
      
  snps = 0;

  // Make a copy (we need to adjust the first base...)

  int errors = 0;

  for (int i=0; i<(int)seq.size(); i++) {
    char c1 = seq[i];
    char c2 = ref[i+startOnRef];
    if (c1 != c2) {
      errors++;
    }
  }

	
  if (bPrint) {
    //cout << "Errors: " << errors << endl;
    cout << "seq ";
    for (int i=0; i<(int)seq.size(); i++) {
      cout << as_base(seq[i]);
    }
    cout << endl << "ref ";
    for (int i=0; i<(int)seq.size(); i++) {
      cout << as_base(ref[i+startOnRef]);
    }
    cout << endl;
  }

  return errors;
}


void CombAligner::UseDoubleSpacing(bool b)
{
  if (b)
    Combmer::UseDoubleSpacing(true);
  else 
    Combmer::UseDoubleSpacing(false);
}



void CombAligner::SetK(int k)
{
  cout << "Setting K=" << k << endl;
  Combmer::SetK((int)k);

}


void CombAligner::MakeCombs()
{
  cout << "Adding combmers..." << endl;
  for (int i=0; i<(int)m_ref.size(); i++) {
    AddCombmers(m_refCombs, m_ref[i], i);
  }
  cout << "Sorting..." << endl;
  Sort(m_refCombs);

}


void CombAligner::ReadFastb(const String & fileName, bool rc)
{
  cout << "Reading reference " << fileName << endl;
  m_ref.ReadAll(fileName);
  cout << "done!" << endl;
  int i, j;
  m_rcStart = (int)m_ref.size();
  if (rc) {
    cout << "rc'ing sequences..." << endl;
    m_ref.resize(m_rcStart * 2);
    for (i=m_rcStart; i<2 * m_rcStart; i++) {
      basevector & b = m_ref[i];
      b = m_ref[i-m_rcStart];
      b.ReverseComplement();
    }
  }
  cout << "Making combs..." << endl;


  int k = Combmer::GetSize();
  
  int total = 0;
  
  for (i=0; i<(int)m_ref.size(); i++) {
    basevector & ref = m_ref[i];
    total += (int)ref.size() - k + 2;
  }
  cout << "Total ref sequences: " << m_ref.size() << endl;
  cout << "Total comb-mers: " << total << endl;
  
  m_refCombs.resize(total);
  int kk = 0;
  for (j=0; j<(int)m_ref.size(); j++) {
    basevector & ref = m_ref[j];
    if (j % 1000 == 0)
      cout << j << endl;
    for (i=0; i<(int)ref.size()-k+2; i++) {
      m_refCombs[kk].Set(ref, i, j, i);
      kk++;
    }  
  }
  cout << "Sorting..." << endl;
  Sort(m_refCombs);
  



  //MakeCombs();
  cout << "done!" << endl;
}

  
void CombAligner::Align(vec<PrelimMatch> & matches, 
			     const basevector & seq,
			     const String & name, 
			     int index)
{
  AlignInt(matches, seq, name, index, true);

}

 
void CombAligner::AlignInt(vec<PrelimMatch> & matches, 
				const basevector & seq,
				const String & name, 
				int index,
				bool bSort)
{
  int j;

  vec<Combmer> seqCombs;
  AddCombmers(seqCombs, seq, 0);
  
  int lastContig = -1;
  int lastStartOnRef = -1;
  
  bool rc;
  if (m_rcStart != -1 && m_rcStart < (int)m_ref.size())
    rc = true;

  for (j=0; j<seqCombs.isize(); j++) {
    
    vec<Combmer>::iterator iter;
    
    vec<Combmer>::iterator begin = m_refCombs.begin();
    vec<Combmer>::iterator end = m_refCombs.end();
    iter = lower_bound(begin, end, seqCombs[j]);
    int pos = distance( begin, iter );
    
    
    while (pos < m_refCombs.isize() && m_refCombs[pos] == seqCombs[j]) {
      int index = m_refCombs[pos].GetStart();
      int indexSeq = seqCombs[j].GetStart();
      
      int startOnRef = index - indexSeq;
      
      //cout << "Testing " << j << " from " << startOnRef << " last=" << lastStartOnRef << endl;
      
      //if (refCombs[pos].GetContig() == lastContig && startOnRef == lastStartOnRef) {
      // pos++;
      //continue;
      //}
      
      lastContig = m_refCombs[pos].GetContig();
      lastStartOnRef = startOnRef;
      
      int errors = 0;
      int snps = 0;
      errors = IsGoodEnough(m_ref[m_refCombs[pos].GetContig()], 
			    seq, 
			    startOnRef,
			    false,
			    snps);
      if (errors <= (int)m_maxErrors) {
	int realRefContig = m_refCombs[pos].GetContig();
	bool rcMatch = false;
	if (rc && realRefContig >= m_rcStart) {
	  realRefContig -= m_rcStart;
	  rcMatch = true;
	}

	if (m_bPrint) {
	  cout << name << " aligns on reference contig " << realRefContig;
	  cout << " at " << startOnRef << " with " << errors << " errors ";
	  if (rcMatch)
	    cout << "rc";
	  else
	    cout << "fw";
	  
	  cout << " and " << snps << " SNPs" << endl;

	  errors = IsGoodEnough(m_ref[m_refCombs[pos].GetContig()], 
				seq, 
				startOnRef,
				true,
				snps);
	  
	  fflush(stdout);
	}
	PrelimMatch tmpM;
	if (rcMatch)
	  startOnRef = (int)m_ref[m_refCombs[pos].GetContig()].size() - startOnRef - (int)seq.size();
	tmpM.Set(index, name, realRefContig, startOnRef, rcMatch, errors, snps);
	matches.push_back(tmpM);	  
	if (m_bPrint) {
	  cout << "---------------------------------------------------" << endl;
	}
      }
      
      
      pos++;
    }
  }

  if (bSort)
    UniqueSort(matches);

}
  




void CombAligner::Align(vec<PrelimMatch> & matches, 
			     const vecbasevector & seq,
			     vec<String> & names)
{

  int i;

  for (i=0; i<(int)seq.size(); i++) {
    AlignInt(matches, seq[i], names[i], i, false);
  }

  UniqueSort(matches);

}
 
