#include "pairwise_aligners/PerfectMatch.h"
#include "Basevector.h"

#include <functional>
#include <set>

 
///Extend a forward perfect match. 
///This method is "private" to the .cc file, and is only meant to be
///called by ExtendPerfectMatch.
void ExtendForwardPerfectMatch(PerfectMatch & m, const basevector & source,
			       const basevector & target);

  /// switch the source and target. We need the lengths if rc1 is true.
void PerfectMatch::Switch(int sourceLength, int targetLength) {
  swap(source, target);
  if (!rc1) {
    swap(startOnSource, startOnTarget);
  } else {
    int temp = startOnTarget;
    startOnTarget = sourceLength - startOnSource;
    startOnSource = targetLength - temp;
  }
}

void PerfectMatch::Canonicalize(int sourceLength, int targetLength) {
  if (source > target ||
      ( source == target && startOnSource > startOnTarget) ) {
    Switch(sourceLength, targetLength);
  }
}

void PerfectMatch::Extend(const basevector & source, 
    const basevector & target) {
  //cout <<"source: ";
  //if (!rc1) source.Print(cout);
  //else { basevector b; b.ReverseComplement(source);b.Print(cout); }
  //cout <<"target: ";
  //target.Print(cout); 
  //PRINT(*this);
  if (rc1) ExtendReversed(source, target);
  else ExtendForward(source, target);
  //PRINT(*this);
}

struct PerfectMatchCompare: 
  public binary_function<bool, const PerfectMatch &,const PerfectMatch &> 
{
  bool operator() (const PerfectMatch & lhs, const PerfectMatch & rhs) const {
    return lhs < rhs;
  }
};

bool operator <(const PerfectMatch & lhs, const PerfectMatch & rhs) {
  if (lhs.source < rhs.source) return true;
  if (lhs.source > rhs.source) return false;
  if (lhs.target < rhs.target) return true;
  if (lhs.target > rhs.target) return false;
  if (lhs.startOnSource < rhs.startOnSource) return true;
  if (lhs.startOnSource > rhs.startOnSource) return false;
  if (lhs.length < rhs.length) return true;
  if (lhs.length > rhs.length) return false;
  if (lhs.startOnTarget < rhs.startOnTarget) return true;
  if (lhs.startOnTarget > rhs.startOnTarget) return false;
  if (lhs.rc1 < rhs.rc1) return true;
  return false;
} 

///Extend a forward perfect match. 
void PerfectMatch::ExtendForward( const basevector & source,
    const basevector & target) {
  Assert(!rc1);
  int old = Offset();
  int backward = 0, forward=0;
  int j = startOnTarget;
  int i = startOnSource;
  for (; i != 0 && j != 0; --i, --j) {
    if (source[i-1] == target[j-1]) --backward;
    else break;
  }

  i = EndOnSource();
  const int MAXI = source.size();

  j = EndOnTarget();
  const int MAXJ = target.size();

  for (; i < MAXI && j < MAXJ; ++i, ++j) {
    if (source[i] == target[j]) ++forward;
    else break;
  }
  startOnSource += backward;
  startOnTarget += backward;
  length += (forward - backward);
  AssertEq(old, Offset());
}

unsigned char reversedBases[4] = {3,2,1,0};
    
///Extend a reversed perfect match. 
void PerfectMatch::ExtendReversedSlow(const basevector & source,
			       const basevector & target) {
  PerfectMatch copy(*this);
  Assert(rc1);
  basevector b;
  b.ReverseComplement(source);
  rc1 = false;
  ExtendForward(b, target);
  rc1 = true;
  copy.ExtendReversed(source, target);
  if (!(copy == *this)) {
    PRINT2(copy, *this);
    FatalErr("slow and fast different");
  }
}

void PerfectMatch::ExtendReversed(const basevector & source,
				      const basevector & target) {
  int backward = 0, forward=0;
  int i = source.size() - startOnSource;
  int j = startOnTarget - 1;
  for (; i < source.isize() && j >= 0; ++i, --j) {
    if (reversedBases[source[i]] == target[j]) --backward;
    else break;
  }

  //cout<<"orig source: "; source.Print(cout);
  i = source.size() - EndOnSource() -1;
  j = EndOnTarget();
  for (; i >= 0 && j < target.isize(); --i, ++j) {
    //PRINT4(i,j,as_base(source[i]), as_base(target[j]));
    if (reversedBases[source[i]] == target[j]) ++forward;
    else break;
  }
  startOnSource += backward;
  startOnTarget += backward;
  length += (forward - backward);
}

    
void ExtendPerfectMatch(PerfectMatch & m, const basevector & source,
                        const basevector & target) {
  if (m.rc1) {
    basevector reverse(source);
    reverse.ReverseComplement();
    m.rc1=false;
    m.ExtendForward(reverse, target);
    m.rc1=true;
  } else {
    m.ExtendForward(source, target);
  }
}

///First, we extend each match
///Then, we remove all matches that are included by any other match.
/// It turns out that we cannot just remove matches that overlap, because
/// of repeats (e.g. CCATATATGG vs. CCATATATATATATGG), and that we also
/// cannot remove only matches that have the same offset, as the previous
/// version did (e.g. CCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATT vs. the same
/// sequence: this will  give many more than one match if we insist on 
/// matching the offsets. However, if one match includes another then since
/// they are both perfect we can be sure that we can remove the second.
///
/// This method is inefficient, but right now I just want to make it work.
/// It checks everything twice and extends everything probably before it
/// needs to.
///
/// This method is used to find alignments. It is not the same method used
/// to find repeats! For repeats, a match included in another should not
/// be absorbed unless the offsets are the same. For alignments on the other
/// hand, a smaller match included in a larger one is likely to be
/// inferior.
void ExtendPerfectMatchesByInclusion(vec<PerfectMatch> & matches,
                          const basevector & source,
                          const basevector & target) {
  int lastIndex = matches.isize()-1;
  for (int i=0; i != matches.isize(); ++i) ExtendPerfectMatch(matches[i],source, target);
  for (int i=0; i <= lastIndex; ++i) {
    for (int j=i+1; j <= lastIndex; ) {
      if (matches[i].Includes(matches[j])) {
      	//PRINT6(i,j,lastIndex,matches[i],"includes",matches[j]);
	swap(matches[j], matches[lastIndex]);
	--lastIndex;
      } else if (matches[j].Includes(matches[i])) {
      	//PRINT6(j,i,lastIndex,matches[j],"includes",matches[i]);
	swap(matches[i], matches[j]);
	swap(matches[j], matches[lastIndex]);
	--lastIndex; 
      } else ++j;
    }
  }
  matches.erase(matches.begin()+lastIndex+1, matches.end());
  //cout << "after inner extension:" << matches << endl;
}
      
