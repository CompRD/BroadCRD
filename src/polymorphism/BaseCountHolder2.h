/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


#ifndef BASE_COUNT_HOLDER2_H
#define BASE_COUNT_HOLDER2_H

#include "CoreTools.h"
#include "Basevector.h"
#include "TokenizeString.h"
#include "math/Functions.h"

#include <numeric>
#include <map>

///Class to keep track of coverage on fw and reverse strands for one contig.
///fw strands are counts_[0..3], while rc are counts_[4..7]
///The class also keeps track of where it is on the contig, via the
///instance variable offsetStart_.
class BaseCountHolder2 {
private:
  typedef map<int, vec<unsigned short> > Cmap;
  typedef Cmap::iterator citer;
  typedef Cmap::const_iterator const_citer;
  Cmap counts_;
  int offsetStart_;
  int contig;

  static const int LINE_SIZE = 8;

public:
  ///Constructor takes fake length argument to be similar to BaseCountHolder.
  BaseCountHolder2(int length=0): counts_(),
				  offsetStart_(0) {}

    
  //void SetNumLines(int maxCount) {}

  int Contig() const { return contig; }
  void SetContig(int c) { contig = c; }

  //const vec<vec<unsigned short> > Counts() const { return counts_; }
  int size() const { return counts_.size(); } 

  void Increment(int offset, unsigned char base, bool rc) {
    vec<unsigned short> * vptr = GetVecAt(offset);
    if (0 == vptr) {
      counts_[offsetStart_+offset] = vec<unsigned short>(8,0);
      //PRINT2(longlong(&(counts_[offsetStart_+offset])), counts_[offsetStart_+offset]);
      vptr = &(counts_[offsetStart_+offset]);
      //PRINT(longlong(vptr));
    }
    //PRINT3(offsetStart_ + offset, contig, (*vptr)[base + rc*4]);
    ++((*vptr)[base + rc*4]);
  }

  void SetOffsetStart( int start) { offsetStart_ = start; }
  int GetOffsetStart() const { return offsetStart_; }
  ///Print one line of data
  void PrintLine(ostream & os , int offset=0) const {
    const  vec<unsigned short> * vptr = GetVecAt(offset);
    if (!vptr) os <<"0\t0\t0\t0\t0\t0\t0\t0\t";
    else copy (vptr->begin(), vptr->end(),
	       ostream_iterator<unsigned short>(os,"\t"));
  }

  int Sum(int offset = 0) const {
    const vec<unsigned short> * vptr = GetVecAt(offset);
    if (!vptr) return 0;
    return accumulate(vptr->begin(), vptr->end(), 0);
  }


  const vec<unsigned short> * GetVecAt(int offset = 0) const {
    offset += offsetStart_;
    const_citer iter = counts_.find(offset);
    if (iter == counts_.end()) return 0;
    else return &(iter->second);
  }

  vec<unsigned short> * GetVecAt(int offset = 0)  {
    offset += offsetStart_;
    citer iter = counts_.find(offset);
    if (iter == counts_.end()) return 0;
    else return &(iter->second);
  }

  pair<int, int> FwRc(int offset =0) {
    offset += offsetStart_;
    citer iter = counts_.find(offset);
    if (iter == counts_.end()) return make_pair(0,0);
    const vec<unsigned short> & d = iter->second;
    return make_pair(d[0]+d[1]+d[2]+d[3],d[4]+d[5]+d[6]+d[7]);
  }


};

#endif //BASE_COUNT_HOLDER2_H
