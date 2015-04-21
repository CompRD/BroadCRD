/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
#include "jumping/ReadFragment.h"
#include "CoreTools.h"

const char * const ReadFragment::typeCodes = "ULGSVC";

std::ostream& operator<<(std::ostream &out, const ReadFragment &rf)
{
  out << "[" << rf.start << " (" << ReadFragment::typeCodes[rf.type] << " "
      << (rf.rc ? "rc " : "fw ")
      << rf.id << " " << rf.target_start << " " << rf.target_end << ") "
      << rf.end << "]";
  return out;
}

bool operator<(const ReadFragment &a, const ReadFragment &b)
{
  if (a.start < b.start) return true;
  if (b.start < a.start) return false;
  if (a.end < b.end) return true;
  if (b.end < a.end) return false;
  if (a.type < b.type) return true;
  if (b.type < a.type) return false;
  if (a.id < b.id) return true;
  if (b.id < a.id) return false;
  if (a.target_start < b.target_start) return true;
  if (b.target_start < a.target_start) return false;
  return (a.target_end < b.target_end);
}

bool operator<(const ReadFragmentHalf &a, const ReadFragmentHalf &b)
{
  if (a.pos < b.pos) return true;
  if (b.pos < a.pos) return false;
  if (a.type < b.type) return true;
  if (b.type < a.type) return false;
  if (a.id < b.id) return true;
  if (b.id < a.id) return false;
  if (a.target_start < b.target_start) return true;
  if (b.target_start < a.target_start) return false;
  return (a.target_end < b.target_end);
}

// A specialized operator== for finding the matching half of a ReadFragmentHalf
bool operator==(const ReadFragmentHalf &a, const ReadFragmentHalf &b)
{
  return (a.seq == b.seq);
}

std::ostream & operator<<(std::ostream &out, const ReadFragmentHalf &rf)
{
  out << "[" << (rf.start ? "" : "<") << rf.pos
      << " (" << ReadFragment::typeCodes[rf.type] << " "
      << rf.id << " " << rf.target_start << " " << rf.target_end << ") "
      << (rf.start ? ">" : "") << "]";
  return out;
}

// TODO: There is an unnecessarily large number of swaps here.  Efficiency could
// be improved greatly. Also, final resize looks just wrong.  TWS.
void RemoveShortFragments(ReadFragmentVec& frag,
			  const int MAX,
			  ReadFragment::TYPE type ) {
  for (ReadFragmentVec::size_type i=0; i < frag.size(); ++i) {
    if (frag[i].type == type && frag[i].Length1() < MAX) {
      for (ReadFragmentVec::size_type j=i; j !=frag.size() -1; ++j) {
	swap(frag[j],frag[j+1]);
      }
      frag.resize(frag.size()-1);
    }
  }
}

#include "feudal/SmallVecDefs.h"
template class SmallVec<ReadFragment,MempoolAllocator<ReadFragment> >;

#include "feudal/OuterVecDefs.h"
template class OuterVec<ReadFragmentVec>;
