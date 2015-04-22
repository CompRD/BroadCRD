#ifndef JUMPING_ANNOTATION_CLEANER_H
#define JUMPING_ANNOTATION_CLEANER_H
/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
#include "jumping/ReadFragment.h"

/// \class AnnotationCleaner Clean up the raw annotations of a read,
/// producing a single annotation for every base on the read (using U
/// annotations where needed).  For use in jumping library analysis to
/// classify the reads based on their content.  
class AnnotationCleaner {
public:
  enum STATE {N=0,///< No annotation
	      L, ///< one L only
	      S, ///< one S only
	      G, ///< one G only
	      GG, ///< multiple G's
	      SG, ///< S and G
	      LG, ///< L and G
	      C}; ///< more than one of S and/or L
  
  vec<ReadFragmentHalf> h; ///< The working vector of all half-ReadFragment's
  vec<ReadFragmentHalf> l, s, g; ///< The stacks of active L,S,G fragments
  int intervalStart; ///<The start of the current annotation interval
  ReadFragmentHalf removed; ///< The most-recently-removed half-ReadFragment
  STATE oldState; ///< The previous state encountered

  AnnotationCleaner() : intervalStart(0) { }
  bool operator()(ReadFragmentVec& annot, int length);
  void markState() { oldState = state(); }
  bool invalid() { return (g.size()>1); }
  void remove(vec<ReadFragmentHalf> &v, const ReadFragmentHalf &r);
  STATE state();
  void emitInterval(ReadFragmentVec& annot, ReadFragmentHalf curr);

};
#endif
