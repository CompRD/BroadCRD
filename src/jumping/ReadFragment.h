#ifndef JUMPING_READ_FRAGMENT_H
#define JUMPING_READ_FRAGMENT_H
/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
#include "system/Types.h"
#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"
#include <ostream>

/// classify a fragment of a read as being of a certain length and type.
/// \class ReadFragment
/// For use in jumping library analysis to classify the reads based on 
/// their content
struct ReadFragment {
  enum TYPE {U=0,//unknown
	     L,//linker 
	     G,//genome
	     S,//stuffer
	     V,//vector
	     C//conflict
  };
  static const char * const typeCodes; ///< One-character codes.  KEEP IN SYNC.
  int start; ///<start of fragment on read, 0-based
  int end; ///<one-past end of fragment on read, 0-based
  TYPE type; ///<what type of alignment, or U if none found
  Bool rc;  ///< whether alignment is forward or rc.
  int id; ///<Which one/which contig of TYPE was found (starting from 0)
  int target_start; ///< Where on the contig the read alignment starts
  int target_end; ///< Where on the contig the read alignment ends

  ReadFragment(int start=0, int end=0, TYPE type=U, Bool rc=False,
	       int id=-1, int target_start=-1, int target_end=-1):
    start(start), end(end), type(type), rc(rc), id(id), target_start(target_start),
    target_end(target_end){}

  int Length1() const { return end - start; }
  int Length2() const { return target_end - target_start; }
};

TRIVIALLY_SERIALIZABLE(ReadFragment);

typedef SerfVec<ReadFragment> ReadFragmentVec;
typedef MasterVec<ReadFragmentVec> VecReadFragmentVec;
extern template class SmallVec<ReadFragment,MempoolAllocator<ReadFragment> >;
extern template class OuterVec<ReadFragmentVec>;

std::ostream& operator<<(std::ostream &out, const ReadFragment &rf);

///compare start, then end, then the other member variables
bool operator<(const ReadFragment &a, const ReadFragment &b);

/// Beginning/ending point of a ReadFragment 
/// \class ReadFragmentHalf
/// For use in jumping library analysis to classify the reads based on 
/// their content TODO: reduce code duplication with ReadFragment
struct ReadFragmentHalf {
  int pos; ///<position of fragment half on read, 0-based
  ReadFragment::TYPE type; ///<what type of alignment, or U if none found
  Bool rc;  ///< whether alignment is forward or rc.
  Bool start; ///< whether this is start or end half of fragment
  int id; ///<Which one/which contig of TYPE was found (starting from 0)
  int target_start; ///< Where on the contig the read alignment starts
  int target_end; ///< Where on the contig the read alignment ends
  int seq; ///<Sequence number used to disambiguate which half-fragment goes with which

  ReadFragmentHalf(const ReadFragment &r = ReadFragment(), bool start=true, int seq=-1):
    pos(start ? r.start : r.end),
    type(r.type), rc(r.rc), start(start),
    id(r.id), target_start(r.target_start), target_end(r.target_end), seq(seq) { }
};

// Sort by position on the read, breaking ties similarly to ReadFragment
bool operator<(const ReadFragmentHalf &a, const ReadFragmentHalf &b);

// A specialized operator== for finding the matching half of a ReadFragmentHalf
bool operator==(const ReadFragmentHalf &a, const ReadFragmentHalf &b);

// Output operator, similar to ReadFragment
std::ostream & operator<<(std::ostream &out, const ReadFragmentHalf &rf);

///Remove unknown fragments shorter than a certain length.
void RemoveShortFragments(ReadFragmentVec& frag,
			  const int MAX,
			  ReadFragment::TYPE type = ReadFragment::U);

#endif
