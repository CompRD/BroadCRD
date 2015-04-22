/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef OVERLAPPING_BAITS_H
#define OVERLAPPING_BAITS_H


#include "Basevector.h"
#include "jumping/ExonStats.h"
#include "Vec.h"

struct Bait {
  explicit Bait(int contig = -1, int begin = -1, int end = -1, int exonId = -1):
    contig(contig), begin(begin), end(end), exonId(exonId) 
  {}

  int contig;
  int begin;
  int end;
  int exonId;

};

inline ostream & operator<<(ostream & os, const Bait & b) {
  os << b.contig << " " << b.begin << " " << b.end << " " << b.exonId << " ";
  return os;
}

inline istream & operator>>(istream & is, Bait & b) {
  is >> b.contig >> b.begin >> b.end >> b.exonId;
  return is;
}

struct SingleBaitParams {
  SingleBaitParams(int fragLength=140):
    fragLength( fragLength)
  { Validate(); }

  int fragLength; 

  ///Assert if invalid parameters.
  void Validate() const;

  int JunkOnEnd() const { return fragLength; }
};

void PickBaits(vec<Bait> & baits, const SingleBaitParams & params,
	       vecbasevector & target, int contig, int start, int length,
	       int exonId = 0);

/// Call PickBaits on each exon and collect all the baits.
void PickAllBaits(vec<Bait> & baits, const SingleBaitParams & params,
		  vecbasevector & ref, const vec<ExonStats> & exons);


#endif // SINGLE_BAITS_H
