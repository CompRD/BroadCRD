/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


#ifndef BASE_COUNT_HOLDER_H
#define BASE_COUNT_HOLDER_H

#include "CoreTools.h"
#include "Basevector.h"
#include "TokenizeString.h"
#include "math/Functions.h"

#include <numeric>

///Class to keep track of coverage on fw and reverse strands for one contig.
///fw strands are counts_[0..3], while rc are counts_[4..7]
///The class also keeps track of where it is on the contig, via the
///instance variable offsetStart_.
class BaseCountHolder {
private:
  vec<vec<unsigned short> > counts_;
  int offsetStart_;
  int contig;

  static const int LINE_SIZE = 8;
public:
  BaseCountHolder(int length = 0): counts_(length, vec<unsigned short>(8,0)),
  offsetStart_(0) {}

  friend istream & operator>>(istream & is, BaseCountHolder b) {
    for (int i = 0; i != LINE_SIZE; ++i) {
      is >> b.counts_[b.offsetStart_];
    }
    return is;
  }
  ///Set one line of data. Note that this method ignores offsetStart_,
  ///because it is meant to be used when reading all the data from a file,
  ///and uses lineno instead. This may be a poor design decision...
  template<class ForwardIterator>
  void SetLine(int lineno, ForwardIterator begin, ForwardIterator end) {
    vec<short unsigned int> & line = counts_[lineno];
    if (line.isize() != LINE_SIZE) line.resize(LINE_SIZE);
    for (int i=0; begin != end; ++begin, ++i) {
      line[i] = *begin;
    }
  }
    
  void SetNumLines(int maxCount) { counts_.resize(maxCount); }

  int Contig() const { return contig; }
  void SetContig(int c) { contig = c; }

  const vec<vec<unsigned short> > Counts() const { return counts_; }
  int size() const { return counts_.size(); } 

  void Increment(int offset, unsigned char base, bool rc) {
    ++counts_[offsetStart_ + offset][base + rc*4];
  }
  unsigned short Get(int offset, unsigned char base, bool rc) const {
    return counts_[offsetStart_ + offset][base + rc*4];
  }
  void SetOffsetStart( int start) { offsetStart_ = start; }
  int GetOffsetStart() const { return offsetStart_; }
  ///Print one line of data
  void PrintLine(ostream & os , int offset=0) const {
    const vec<unsigned short> & data = counts_[offsetStart_ + offset];
    copy (data.begin(), data.end(),ostream_iterator<unsigned short>(os,"\t"));
  }
  int Sum(int offset = 0) const {
    const vec<unsigned short> & data = counts_[offsetStart_ + offset];
    return accumulate(data.begin(), data.end(), 0);
  }
  longlong Total() const {
    longlong ret=0;
    for (int i=0; i != counts_.isize(); ++i) {
      ret += ::Sum(counts_[i]);
    }
    return ret;
  }
  int Uncovered() const {
    int ret=0;
    for (int i=0; i != counts_.isize(); ++i) {
      if (0 == ::Sum(counts_[i])) ++ret;
    }
    return ret;
  }
  int OneCoveredOrLess() const {
    int ret=0;
    for (int i=0; i != counts_.isize(); ++i) {
      if (1 >= ::Sum(counts_[i])) ++ret;
    }
    return ret;
  }
  pair<int, int> FwRc(int offset =0) const {
    const vec<unsigned short> & d = counts_[offsetStart_ + offset];
    return make_pair(d[0]+d[1]+d[2]+d[3],d[4]+d[5]+d[6]+d[7]);
  }
  pair<int, int> Alleles(int offset =0) const {
    const vec<unsigned short> & d = counts_[offsetStart_ + offset];
    vec<int> al(4,0);
    for (int i=0; i !=4; ++i) al[i]=d[i]+d[i+4];
    sort(al.begin(), al.end());
    return make_pair(al[2],al[3]);
  }
};
///Load in the section of a coverage file corresponding to one contig.
///Put the raw coverage in raw, the accepted coverage in accepted,
/// and the reference sequence in ref if it is not NULL.
/// return the contig number or -1 if the istream has no data.
int LoadCoverageContig(istream & is, BaseCountHolder & raw,
		       BaseCountHolder & accepted, basevector * ref = NULL);

///Load in a coverage file.
///Put the raw coverage in raw, the accepted coverage in accepted,
/// and the reference sequence in ref if it is not NULL.
/// return the number of contigs
/// Calls LoadCoverageContig
int LoadCoverageFile(const String & fname, vec<BaseCountHolder> & raw,
		      vec<BaseCountHolder> & accepted, 
		      vecbasevector * ref = NULL);


#endif //BASE_COUNT_HOLDER_H
