// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

/** Remove all unique elements in a sorted container.
    Works in linear time and with given space.
\file RemoveUnique.h
*/

#ifndef REMOVE_UNIQUE_H
#define REMOVE_UNIQUE_H

//Not defining NDEBUG because the Assert(is_sorted) should stay unless 
//we really do not need it.

#include "STLExtensions.h"

template<class T> void MarkAsBad(T & t) { t.MarkAsBad(); }
template<class T> bool IsBad(T & t) { return t.IsBad(); }

/// Equivalent to remove_if(begin,end, is_different_from_previous_and_next);
/// Linear time and operates in place: no extra space needed.
///
/// Caution is needed if working with a set or map, 
/// because this function modifies elements in place.
///
/// Prerequisite: the range (begin, end) is sorted with <, and 
/// operator< sorts into equivalence classes according to Equal.
template<class ForwardIter, class Equal>
ForwardIter RemoveUnique(const ForwardIter & begin, const ForwardIter & end, 
		       Equal eq) {

  if (begin==end) return end;
  Assert(is_sorted(begin, end));
  ForwardIter prev = begin, current = begin, next = begin;
  ++current;
  if (current == end) return begin;//only one element, so remove it!
  ++(++next);

  //deal with the first element; we know there are at least two elements.
  if (! eq(*current,*prev) ) MarkAsBad(*prev);
  //the general case
  while (next != end) {
    if (! (eq(*current, *prev) || eq(*current, *next)) ) MarkAsBad(*current);
    ++prev; ++current; ++next;
  }
  //now deal with the last element.
  if (! eq(*current, *prev) ) MarkAsBad(*current);  
  return remove_if(begin, end, IsBad<typename ForwardIter::value_type>);
}

#endif //REMOVE_UNIQUE_H
