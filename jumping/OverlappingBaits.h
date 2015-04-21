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

struct BaitPair {
  explicit BaitPair(int contig = -1, int fw = -1, int rc = -1,
		    int over = -1, float tm = -1, float gcPct = -1,
		    int homopol = -1, int exonId = -1):
    contig(contig), fwStart(fw), rcStart(rc), overlap(over), tm(tm),
    gcPct(gcPct), homopol(homopol), exonId(exonId) 
  {}

  int contig;
  int fwStart;
  int rcStart;
  int overlap;
  float tm, gcPct;
  int homopol;
  int exonId;

  friend bool operator<(const BaitPair& l,  const BaitPair& r) { 
    return l.tm < r.tm; 
  } 

  ///This is not always equal to fwStart + fragmentlength.
  int fwLength() { return rcStart + overlap - fwStart; }
};

inline ostream & operator<<(ostream & os, const BaitPair & b) {
  os << b.contig << " " << b.fwStart << " " << b.rcStart << " " 
     << b.overlap << " " << b.tm << " " << b.gcPct << " "
     << b.homopol << " " << b.exonId << " ";
  return os;
}

inline istream & operator>>(istream & is, BaitPair & b) {
  is >> b.contig >> b.fwStart >> b.rcStart >> b.overlap >> b.tm >> b.gcPct
     >> b.homopol >> b.exonId;
  return is;
}

struct BaitParams {
  BaitParams(float tm=60.0, float tmMin = 55.0, float tmMax = 65.0,
	     float gcMin = 20.0, float gcMax = 80.0, int homopolMax = 5,
	     float NaMol = 0.05, int fragLength=110, int startWiggle=5, 
	     int minOverlap=20,int maxOverlap=45):
    tm(tm), tmMin(tmMin), tmMax(tmMax), gcMin(gcMin), gcMax(gcMax),
    homopolMax(homopolMax),
    NaMol(NaMol),
    fragLength( fragLength),
    startWiggle(startWiggle),
    minOverlap( minOverlap),
    maxOverlap(maxOverlap)
  { Validate(); }

  float tm, tmMin, tmMax, gcMin, gcMax;
  int homopolMax;
  float NaMol;
  int fragLength; 
  int startWiggle;
  int minOverlap; 
  int maxOverlap;
  int JunkOnEnd() const { return 2 * fragLength  - minOverlap; }

  ///Assert if invalid parameters.
  void Validate() const;
};

struct BaitPairGood {
  const BaitParams &pars;
  BaitPairGood(const BaitParams &pars) : pars(pars) { }
  bool operator()(const BaitPair &p)
  { return (pars.tmMin <= p.tm  && p.tm <= pars.tmMax)
      && (pars.gcMin <= p.gcPct && p.gcPct <= pars.gcMax)
      && (p.homopol <= pars.homopolMax); }
};


/// Given a target basevector, pick baits in pairs that will
/// cross-hybridize. We assume that the target includes, at the end
/// of it, some extra sequence of length at least 
/// params.JunkOnEnd() that is junk sequence: that is,
/// it can be included in a bait if needed, but does not need to
/// be included.
/// We will pick baits with overlaps as 
/// close to params.tm as possible.
/// The output is in baits. reverse baits are
/// always of length params.fragLength, forwardBaits can be up to 
/// params.startWiggle shorter if it's better to start before the exon.
/// The output vector is not cleared, and is added to with 
/// push_back. This allows one to call this function repeatedly 
/// with the same output vector and different targets (e.g. different exons).
//// Returns the number of baits picked for this target.
int PickBaits(vec<BaitPair> & baits, const BaitParams & params,
	       vecbasevector & target, int contig, int start, int length,
	       int exonId = 0);

/// Call PickBaits on each exon and collect all the baits.
/// We lengthen each exon by params.JunkOnEnd() on the right
/// as required by PickBaits.
/// Returns the number of targets for which a bait was picked. 
int PickAllBaits(vec<BaitPair> & baits, const BaitParams & params,
		  vecbasevector & ref, const vec<ExonStats> & exons);


float gcPercent(basevector &b);
int maxHomopolymerLength(basevector &b);

#endif // OVERLAPPING_BAITS_H
