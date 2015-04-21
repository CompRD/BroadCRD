// Copyright (c) 2004-5 Broad Institute/Massachusetts Institute of Technology

#ifndef PERFECT_MATCH_H
#define PERFECT_MATCH_H


#include "system/Types.h"
#include "Basevector.h"
#include "RemoveUnique.h"
#include "STLExtensions.h"
#include "Vec.h"

#include <iostream>

class PerfectMatch;
inline ostream & operator<<(ostream & os, const PerfectMatch & m);

///Represents a perfect match between two sequences.
/// \class PerfectMatch
///
///If rc1 is true, then startOnSource refers to 
///the start position on the ReverseComplement(), that is counting from
///the back of the original source.
///
/// Note that this class is not a POD (plain old data-type) according to
/// the standard because it has a constructor, but in fact it should
/// behave like a POD and be bit-copyable.
///
/// for a PerfectMatch to be rc1 always means that the 
/// source has been reversed, and that
/// we count backwards from the end of the source basevector.
struct PerfectMatch {
  PerfectMatch(int s=0, int ss=0, int t=0, int st=0, int l=0, Bool r=false):
    source(s), startOnSource(ss), target(t), startOnTarget(st),
    length(l), rc1(r) {}
  int source; ///<identifier for the source (e.g. read id).
  int startOnSource; 
  int target; ///<identifier for the target (e.g. contig id).
  int startOnTarget;
  int length; ///<length of the perfect match.
  Bool rc1; ///<whether source is reversed in order to match target.
  ///If so, then it is important to remember that startOnSource refers to 
  ///the start position on the ReverseComplement(), that is counting from
  ///the back of the original source.

  ///difference between the starts.
  int Offset() const { 
    return startOnTarget - startOnSource;
  }
  int EndOnTarget() const { return startOnTarget + length; }
  int EndOnSource() const { return startOnSource + length; }

  /// switch the source and target. We need the lengths if rc1 is true.
  void Switch(int sourceLength, int targetLength);

  /// ensure that source <= target and startOnSource <= startOnTarget.
  void Canonicalize(int sourceLength, int targetLength);

  ///Calls Canonicalize(int, int)
  void Canonicalize(const vec<int> & sourceLengths, 
		    const vec<int> & targetLengths) {
    Canonicalize(sourceLengths[source], targetLengths[target]);
  }

  ///Prerequisite: offsets are equal.
  ///Prerequisite: *this has been extended.
  bool OverlapsThisExtended(const PerfectMatch & o) const {
    AssertEq(Offset(), o.Offset());
    return (o.startOnSource >= startOnSource &&
	    o.startOnSource <= EndOnSource() );
  }
      
  ///Extend a forward match as far as possible, forward and backward.
  ///This should conserve Offset(): should assert that.
  void ExtendForward(const basevector & source, const basevector & target);

  ///Extend a rc1 match as far as possible, forward and backward.
  ///This should conserve Offset(): should assert that.
  void ExtendReversed(const basevector & source, const basevector & target);

  ///Extend a rc1 match as far as possible, forward and backward.
  ///This should conserve Offset(): should assert that.
  void ExtendReversedSlow(const basevector & source, const basevector & target);

  ///Extend as far as possible forward and back, considering our orientation.
  ///This should conserve Offset(): should assert that.
  void Extend(const basevector & source, const basevector & target);

  ///Extend as far as possible forward and back, considering our orientation.
  ///Convenience function, calls Extend(basevector, basevector).
  void Extend(const vecbasevector & source, const vecbasevector & target) {
    Extend(source[this->source], target[this->target]);
  }

  ///Set source to -1
  void MarkAsBad() { source = -1; }

  ///Check if source is -1
  bool IsBad() const { return (-1 == source); }

  ///True if they overlap on both source and target.
  bool Includes (const PerfectMatch & o) const {
    return (source == o.source && 
	    target == o.target &&
	    startOnSource <= o.startOnSource &&
	    startOnTarget <= o.startOnTarget &&
	    EndOnSource() >= o.EndOnSource() &&
	    EndOnTarget() >= o.EndOnTarget());
  }

  ///Like the name says... this is a good order when extending/compacting.
  bool LessByIdOffsetReversedStart(const PerfectMatch & o) const {
    if (source < o.source) return true; 
    if (source > o.source) return false; 
    if (target < o.target) return true; 
    if (target > o.target) return false; 
    if (Offset() <  o.Offset()) return true;
    if (Offset() >  o.Offset()) return false;
    if (rc1 < o.rc1) return true; 
    if (rc1 > o.rc1) return false; 
    if (startOnSource < o.startOnSource) return true;
    return false;
  }

  /// If this is true, then overlap/compaction is possible, but
  /// still requires OverlapsThisExtended().
  bool EqualByIdOffsetReversed(const PerfectMatch & o) const {
    if (source != o.source) return false; 
    if (target != o.target) return false; 
    if (Offset() != o.Offset()) return false;
    if (rc1 != o.rc1) return false; 
    return true;
  }

};


struct SortPerfectMatchByOffset: 
  public binary_function<PerfectMatch, PerfectMatch, bool> {
  bool operator() (const PerfectMatch & l, const PerfectMatch & r) {
    return l.Offset() < r.Offset();
  }
};

struct SortPerfectMatchForExtending:
  public binary_function<PerfectMatch, PerfectMatch, bool> {
  bool operator() (const PerfectMatch & lhs, const PerfectMatch & rhs) {
    return lhs.LessByIdOffsetReversedStart(rhs);
  }
};

inline bool operator==(const PerfectMatch & rhs, const PerfectMatch & lhs) {
  return (rhs.source == lhs.source && rhs.startOnSource == lhs.startOnSource
      && rhs.target == lhs.target && rhs.startOnTarget == lhs.startOnTarget
      && rhs.length == lhs.length && rhs.rc1 == lhs.rc1);
}

inline ostream & operator<<(ostream & os, const PerfectMatch & m) {
  os << m.source << " " << m.startOnSource << " " << m.target << " " 
    << m.startOnTarget << " " << m.length << " " 
     << (m.rc1 ? "rc1" : "fw");
  return os;
}

///Sort by source, then by target, then by start on target, then by length.
bool operator<(const PerfectMatch & rhs, const PerfectMatch & lhs);

///Extend one perfect match as far as it can go.
///We assume that the match actually belongs to the given source and target!
void ExtendPerfectMatch(PerfectMatch & m, const basevector & source,
                        const basevector & target);

///Extend a vector of perfect matches as far as they can go.
///Note that the vector itself is modified, so if you pass in 20 matches
/// that can all be aligned together, you will get back a vector with 
/// only one match.
/// The method tries to be efficient, and checks whether a match is
/// included in a prior extension before trying to extend it.
///
/// This method will remove matches if they are included in a
/// previously extended match even if they have a different offset. Note that
/// this differs from ExtendAndRemoveOverlaps.
///
/// This is meant to be used for things like FixLocationsWithKmers
void ExtendPerfectMatchesByInclusion(vec<PerfectMatch> & m, 
				     const basevector & source,
				     const basevector & target);


/// Extend a vector of perfect matches as far as they can go.
///
/// This method works like (and uses) remove_if: you get back an iterator that
/// points one-past the last "good" match, and you can then happily 
/// erase from that return point to end.
///
/// This method will only remove matches if they are included in a
/// previously extended match and also have the same offset. Note that
/// this differes from ExtendPerfectMatchesByInclusion.
///
/// This is meant to be used for things like TagRepeats, where we want
/// to find small repeats even if they are included inside larger ones.
template<class ForwardIter>
ForwardIter ExtendAndRemoveOverlaps(ForwardIter begin, ForwardIter end,
				    const vecbasevector & source, 
				    const vecbasevector & target) {
  if (begin == end) return end;

  SortPerfectMatchForExtending sorter; 
  if (!is_sorted(begin, end, sorter)) {
    sort(begin, end, sorter);
  }

  ForwardIter previous = begin, current = begin;
  previous->Extend(source, target);

  for (++current; current != end; ++current ) {
    if (previous->EqualByIdOffsetReversed(*current)
	&& previous->OverlapsThisExtended(*current) ) {
      current->MarkAsBad();
    }
    else {
      current->Extend(source, target);
      previous = current;
    }
  }
  return remove_if(begin, end, IsBad<typename ForwardIter::value_type>);
}  
  

#endif //PERFECT_MATCH_H
