/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
#include "CoreTools.h"
#include "jumping/AnnotationCleaner.h"

bool AnnotationCleaner::operator()(ReadFragmentVec& annot, int length)
{
  ReadFragmentVec out;
  // Make sure nothing is hanging around
  h.clear();
  l.clear();
  s.clear();
  g.clear();
  intervalStart = 0;
  // Transform ReadFragment's into ReadFragmentHalf's
  // All the needed information is now in h; we will completely rewrite annot
  int seq = 0; // sequence number used to disambiguate pairing of halves
  for (ReadFragmentVec::size_type i=0; i<annot.size(); ++i) {
    h.push_back(ReadFragmentHalf(annot[i], true, seq));
    h.push_back(ReadFragmentHalf(annot[i], false, seq));
    ++seq;
  }
  // Put the fragments in order along the sequence.
  stable_sort(h.begin(), h.end());
  // Now walk through, keeping stacks of the fragments we're in the
  // middle of.  The genome classification is weakest, yielding to
  // stuffer which yields to linker.
  for (vec<ReadFragmentHalf>::iterator curr = h.begin(); curr != h.end() ; ++curr) {
    markState();
    if (curr->start) {
      // Record interval appropriately
      if (curr->type==ReadFragment::L) 
	l.push_back(*curr);
      else if (curr->type==ReadFragment::S)
	s.push_back(*curr);
      else if (curr->type==ReadFragment::G)
	g.push_back(*curr);
      emitInterval(out, *curr);
    } else {
      // Remove appropriate interval
      if (curr->type==ReadFragment::L) 
	remove(l, *curr);
      else if (curr->type==ReadFragment::S)
	remove(s, *curr);
      else if (curr->type==ReadFragment::G)
	remove(g, *curr);
      emitInterval(out, *curr);
    }
    if (invalid()) {
      PRINT((*curr));
      // Invalid read (repetitive) detected!  Punt.
      annot.clear();
      // emit an ambiguous G interval covering the entire read
      annot.push_back(ReadFragment(0, length, ReadFragment::G, false));
      return false; // RETURN!
    }
  }
  // One more interval may be needed to cover from last alignment to end of read
  markState();
  emitInterval(out, ReadFragmentHalf(ReadFragment(length, length)));
  annot = out;
  return true;
}

void AnnotationCleaner::remove(vec<ReadFragmentHalf> &v, const ReadFragmentHalf &r)
{
  vec<ReadFragmentHalf>::iterator it = std::find(v.begin(), v.end(), r);
  ForceAssert(it!=v.end());
  removed = *it;
  it = std::remove(it, v.end(), r);
  ForceAssertEq(1, distance(it, v.end()));
  v.erase(it, v.end());
}

AnnotationCleaner::STATE AnnotationCleaner::state()
{
  unsigned int nannot = l.size() + s.size() + g.size();
  if (nannot==0)
    return N;
  if (nannot==1) {
    if (l.size()>0)
      return L;
    if (s.size()>0)
      return S;
    return G;
  }
  // nannot >= 2
  if (g.size()==nannot)
    return GG;
  if (g.size()+1==nannot) {
    if (s.size()==1)
      return SG;
    if (l.size()==1)
      return LG;
  }
  return C;
}

void AnnotationCleaner::emitInterval(ReadFragmentVec& annot,
				     ReadFragmentHalf curr)
{
  STATE newState = state();
  // The interval's type will get filled in by cases below.  Typically
  // when we emit, curr is a close, and its id/start/end are what we
  // want; occasionally we emit based on an open, and then need to go
  // find the right interval to emit.
  ReadFragment interval(intervalStart, curr.pos, ReadFragment::U,
			curr.rc, curr.id, curr.target_start, curr.target_end); 
  if ((oldState==S && (newState==C || newState==N))
      || (oldState==SG && newState==G)) {
    interval.type = ReadFragment::S;
    if (newState==C) {
      // The previous S state, not curr, is what we need for our interval
      ForceAssertGe(s.size(), 1u);
      interval.id = s[0].id;
      interval.target_start = s[0].target_start;
      interval.target_end = s[0].target_end;
    }
  } else if ((oldState==L && (newState==C || newState==N))
	     || (oldState==LG && newState==G)) {
    interval.type = ReadFragment::L;
    if (newState==C) {
      // The previous L state, not curr, is what we need for our interval
      ForceAssertGe(l.size(), 1u);
      interval.id = l[0].id;
      interval.target_start = l[0].target_start;
      interval.target_end = l[0].target_end;
    }
  } else if (oldState==G) {
    ForceAssert(newState==N || newState==SG || newState==LG || newState==GG);
    interval.type = ReadFragment::G;
    if (newState != N) {
      ForceAssertGe(g.size(), 1u);
      interval.id = g[0].id;
      interval.target_start = g[0].target_start;
      interval.target_end = g[0].target_end;
    }
  } else if (oldState==C && (newState==S || newState==L)) {
    interval.type = ReadFragment::C;
    interval.id = -1;
    interval.target_start = -1;
    interval.target_end = -1;
  } else if (oldState==N) {
    ForceAssert(newState==S || newState==L || newState==G
		|| oldState==N); // the last call is an N-to-N transition.
    interval.type = ReadFragment::U;
    interval.id = -1;
    interval.target_start = -1;
    interval.target_end = -1;
  } else
    return; // RETURN!
  if (interval.start != interval.end) {
    annot.push_back(interval);
    intervalStart = interval.end; // update the start position of next interval
  }
  return;
}

