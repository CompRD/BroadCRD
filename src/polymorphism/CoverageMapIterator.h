#ifndef COVERAGE_MAP_ITERATOR_H
#define COVERAGE_MAP_ITERATOR_H
/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// Class to provide a one-line-at-a-time view of a coverage map file.

#include "CoreTools.h"
#include "TokenizeString.h"
#include "polymorphism/BaseMapInfo.h"

class ContigCoverageIterator;

class CoverageMapIterator {
  ifstream in;
  vec<String> x;
  string line;
  int contig_;
  int nBases_;
  bool done_;
  BaseMapInfo bases_;
  /// Initialize contig and nbases information, and skip over header line
  void InitializeContig()
  {
    Tokenize(line, x, " ");
    if (in && line.size()>0 && x[0]=="Contig") {
      contig_ = x[1].Int();
      nBases_ = x[3].Before(":").Int();
      getline(in, line); // Skip header
      if (!in)
	done_ = true;
    } else {
      done_ = true;
    }
  }

public:
  
  CoverageMapIterator(String filename) : in(filename.c_str()), done_(false)
  {
    getline(in, line);
    InitializeContig();
    this->operator++();
  }
  
  /// Move to next data line, updating contig if we hit a new one.
  void operator++()
  {
    while (!done_) {
      getline(in, line);
      if (!in)
	done_ = true;
      else {
	x.clear();
	Tokenize(line, x, "\t");
	if (x.size() < 18) {
	  // This is probably a new Contig line
	  InitializeContig();
	} else {
	  bases_ = BaseMapInfo(x);
	  break; // We actually found another data line...
	}
      }
    }
  }
  
  const BaseMapInfo & operator*() { return bases_; }
  const BaseMapInfo * operator->() { return &bases_; }

  friend bool done(const CoverageMapIterator &it) { return it.done_; }
  friend int nBases(const CoverageMapIterator &it) { return it.nBases_; }
  friend int contig(const CoverageMapIterator &it) { return it.contig_; }

};

class ContigCoverageIterator {
  CoverageMapIterator &it_;
  const int contig_;
  bool done_;
public:
  ContigCoverageIterator(CoverageMapIterator &it) : it_(it),
						    contig_(contig(it)),
						    done_(false)
  { }
  void operator++()
  {
    ++it_;
    if (done(it_) || contig(it_) != contig_)
      done_=true;
  }
  const BaseMapInfo & operator*() { return *it_; }
  const BaseMapInfo * operator->() { return (it_.operator->()); }
  friend bool done(const ContigCoverageIterator &it)  { return it.done_; }
  
};

inline ContigCoverageIterator contigIt(CoverageMapIterator &it)
{
  return ContigCoverageIterator(it); 
}


#endif //COVERAGE_MAP_ITERATOR_H
