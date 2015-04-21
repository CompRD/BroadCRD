///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LOOK_ALIGN_FINDER_H
#define LOOK_ALIGN_FINDER_H

#include "lookup/LookAlign.h"
#include "Vec.h"

/**
   Class: LookAlignFinder

   Iterate through a <.qltout> file, grabbing all aligns for one read each time.

   The point is to use as little memory as possible but still get all the
   alignments for each read.

   We assert that the .qltout file is not necessarily ordered, but is clumped:
   that is, all aligns for one read are consequtive.
*/
class LookAlignFinder {
private:
  vec<look_align_plus> aligns_;
  look_align_plus last_;
  set<int> found_;
  ifstream * isp_;

  void SetEmpty() { last_.query_id = -1; }

public:

  /// Constructor: LookAlignFinder constructor
  /// Initialize from a .qltout file.
  LookAlignFinder(const String & fname);

  ///Close the input stream and clean up.
  ~LookAlignFinder() { isp_->close(); delete isp_; }

  /// Method: QueryId
  /// Read number for the current set of alignments.
  int QueryId() const { return aligns_.empty() ? -1 : aligns_[0].query_id; }


  /// Method: Aligns constant
  /// Current vector of alignments (not sorted).
  const vec<look_align_plus> & Aligns() const { return aligns_; }

  /// Method: Aligns
  /// Current vector of alignments (not sorted).
  /// *Be very careful* with this function! You can sort them to your
  /// heart's content, but if you change the query_ids you will 
  /// create serious errors in operator++ and probably lead to a crash.
  vec<look_align_plus> & Aligns()  { return aligns_; }

  LookAlignFinder & operator++();

  bool empty() const { return -1 == last_.query_id; }
};



#endif //LOOK_ALIGN_FINDER_H
