/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef COMB_ALIGNER_
#define COMB_ALIGNER_

#include "Basevector.h"


#include "String.h"
#include "Vec.h"
#include <ostream>




class Combmer
{
public:
  Combmer() {
    m_contig = -1;
    m_start = -1;
    //m_k = 10;
    //m_skip = 2;
    //m_bDouble = true;
  }

  ~Combmer() {}


  void Print() const {
    cout << "Combmer, c" << m_contig << " -> ";
    for (int i=0; i<(int)m_bases.size(); i++) {
      cout << (int)m_bases[i];
    }
    cout << endl;
  }

  static int GetSize() {
    if (Combmer::m_bDouble == true) {
      return Combmer::m_k * Combmer::m_skip - 1;      
    } else {
      return Combmer::m_k * Combmer::m_skip;      
    }
  }
  static void SetK(int k) {
    Combmer::m_k = k;
    cout << "Using K=" << Combmer::m_k << endl;
  }

  static void UseDoubleSpacing(bool b)
  {
    Combmer::m_bDouble = b;
  }

  bool Set(const basevector & c, 
	   int offset,
	   int contig,
	   int start) {

 
    m_start = start;
    m_contig = contig;

    int i;
    m_bases.resize(Combmer::m_k);
    

    if (m_bDouble) {
      int plus = 0;
      for (i=0; i<Combmer::m_k; i++) {
	int pos = offset + i + plus;
	char col = c[pos];
	m_bases.Set(i, col);	
	if (i % 2 == 1)
	  plus += 2;
      }
    } else {
      for (i=0; i<Combmer::m_k; i++) {
	int pos = offset + Combmer::m_skip * i;
	if (pos >= (int)c.size())
	  break;
	char col = c[pos];
	m_bases.Set(i, col);	
      }
    }
    return true;
  }

  int GetContig() const {return m_contig;}
  int GetStart() const {return m_start;}

  bool operator < (const Combmer & c) const {
    if (m_bases == c.m_bases) {
      if (m_contig != c.m_contig)
	return (m_contig < c.m_contig);
      return (m_start < c.m_start); 
    }
    return m_bases < c.m_bases;
  }
  bool operator == (const Combmer & c) const {
    return (m_bases == c.m_bases);
  }

private:
  basevector m_bases;
  int m_contig;
  int m_start;
  static bool m_bDouble;

  static int m_k;
  static int m_skip;
};




class PrelimMatch
{
public:
  PrelimMatch() {
    m_index;
    m_contig;
    m_start;
    m_rc;
    m_errors = 0;
    m_snps = 0;
  }

  void Set(int index, const String & name, int contig, int start, bool rc, int errors, int snps) {
    m_name = name;
    m_index = index;
    m_contig = contig;
    m_start = start;
    m_rc = rc;
    m_errors = errors;
    m_snps = snps;
  }

  bool operator < (const PrelimMatch & c) const {
    if (m_name != c.m_name) {
      return m_name < c.m_name;
    }
    if (m_contig != c.m_contig) {
      return (m_contig < c.m_contig);
    }
    if (m_rc != c.m_rc) {
      if (m_rc)
	return false;
      else
	return true;
    }
    return (m_start < c.m_start);
  }

  bool operator == (const PrelimMatch & c) const {
    if (m_name != c.m_name) {
      return false;
    }
    if (m_contig != c.m_contig) {
      return false;
    }
    if (m_rc != c.m_rc) {
      return false;
    }
    if (m_start != c.m_start)
      return false;
    return true;
  }

  friend std::ostream& operator<<( std::ostream& out, PrelimMatch const& pm ) {
    out << pm.m_name << '\t' << pm.m_contig << '\t' << pm.m_start << '\t'
            << (pm.m_rc ? '+' : '-') << '\t' << pm.m_errors << '\t'
            << pm.m_snps << '\n';
    return out;
  }
 
  const String & GetName() const {return m_name;}
  int GetRefContig() const {return m_contig;}
  int GetRefStart() const {return m_start;}
  bool IsRC() const {return m_rc;}
  
  int GetNumErrors() const {return m_errors;}

  // CURRENTLY NOT IN USE:
  int GetNumSNPs() const {return m_snps;}

private:
  String m_name;
  int m_index;
  int m_contig;
  int m_start;
  bool m_rc;
  int m_errors;
  int m_snps;
};




//--------------------------------------------------------------------

class CombAligner
{
 public:
  CombAligner() {
    m_maxErrors = 4;
    m_bPrint = true;
    m_bSort = true;
    m_rcStart = -1;
  }

  ~CombAligner() {}


  void UseDoubleSpacing(bool b);

  void SetK(int k = 10);

  void SetPrintAligns(bool b = true) {
    m_bPrint = b;
  }

  void SetMaxErrors(int max) {
    m_maxErrors = max;
  }

  // Reads and converts the reference. If rc == true, it will also
  // add the reverse complement, otherwise only fw matches will
  // be returned.
  void ReadFastb(const String & fileName, bool rc = true);

  // Access to the reference object, so you can read it yourself
  // if you so desire.
  vecbasevector & Reference() {return m_ref;}
  void MakeCombs();


  void Align(vec<PrelimMatch> & matches, 
	     const basevector & seq, 
	     const String & name = "sequence", 
	     int index = 0);

 
  void Align(vec<PrelimMatch> & matches, 
	     const vecbasevector & seq,
	     vec<String> & names);


 private:
  void AlignInt(vec<PrelimMatch> & matches, 
		const basevector & seq, 
		const String & name, 
		int index,
		bool bSort);

  vecbasevector m_ref;
  
  int m_maxErrors;
  int m_rcStart;
  bool m_bPrint;
 
  vec<Combmer> m_refCombs;
  bool m_bSort;
};




#endif //COMB_ALIGNER_

