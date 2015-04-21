/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "dna/DNAHybridization.h"
#include "jumping/OverlappingBaits.h"
#include "math/Functions.h"



void BaitParams::Validate() const {
  ForceAssertLt(maxOverlap, fragLength);
  ForceAssertLt(minOverlap, maxOverlap);
  ForceAssertLt(2*startWiggle, minOverlap);//overlaps between baits must
  //be less than the hybridization overlap.
}

///Modify this struct if you want to use a different picking method.
struct BetterPair: public binary_function<bool, BaitPair, BaitPair> {
  const double target;
  BetterPair(double tm): target(tm) {}
  bool operator()( const BaitPair & l, const BaitPair & r) {
    return abs(l.tm - target) < abs(r.tm - target);
  }
};


float gcPercent( basevector &b )
{
    unsigned int gc = 0;
    unsigned int nnn = b.size();
    for ( unsigned int i = 0; i < nnn; ++i )
        if ( IsGC(b[i]) )
            gc += 1;
    return 100.f * gc / nnn;
}

int maxHomopolymerLength(basevector &b) {
  unsigned char curr = b[0];
  int run=1, max_run=0;
  for (int i=1; i < int(b.size()); ++i) {
    if (curr==b[i]) {
      ++run;
    } else {
      curr = b[i];
      if (run>max_run)
	max_run=run;
      run=1;
    }
  }
  return max(run, max_run);
}

int PickBaits(vec<BaitPair> & baits, const BaitParams & pars,
	       vecbasevector & target, int contig,
	       int targetStart, int targetLength, int exonId ) {
  pars.Validate();

  ForceAssertGt(targetLength, pars.JunkOnEnd());

  int count = 0;
  int baitstart = max(0,targetStart - pars.startWiggle);
  int lastBaitStart = targetStart;
  BetterPair comp(pars.tm);//to pick pair with closest Tm.
  BaitPairGood ok(pars);

  while (baitstart < targetStart + targetLength - pars.JunkOnEnd()) {
    vec<BaitPair> pairs;
    //try all forward fragments within startWiggle.
    for (int start = baitstart; start != lastBaitStart; ++start) {
      int end = start + pars.fragLength;
      //try all rc fragments with acceptable overlaps.
      for (int j = end - pars.maxOverlap; j != end - pars.minOverlap; ++j) {
	basevector overlap;
	overlap.SetToSubOf(target[contig], j, end - j);
	BaitPair site(contig, start, j,
		      overlap.size(),
		      Tm(overlap, pars.NaMol),
		      gcPercent(overlap),
		      maxHomopolymerLength(overlap),
		      exonId);
	if (ok(site)) {
	  pairs.push_back(site);
	}
      }
    }
    if (!pairs.empty()) {
      baits.push_back(*min_element(pairs.begin(), pairs.end(), comp));
      ++count;

      //ensure an overlap between bait pairs of [1,1 + startWiggle);
      baitstart = baits.back().rcStart + pars.fragLength
	- (1 + pars.startWiggle);
      //don't run off the end of the contig.
      lastBaitStart = min(baitstart +  pars.startWiggle,
			  target[contig].isize() - pars.JunkOnEnd());
    } else {
      cout << "Bait design failed in range " << baitstart << " to "
	   << lastBaitStart << " for exon " << exonId << " on contig "
	   << contig << " at " << targetStart << " length " << targetLength
	   << "\n";
      // Bump over and try again
      baitstart = lastBaitStart;
      lastBaitStart += pars.startWiggle;
    }
  }
  return count;
}

int PickAllBaits(vec<BaitPair> & baits, const BaitParams & pars,
		  vecbasevector & target, const vec<ExonStats> & exons) {
  baits.clear();
  int count=0;
  for (int i=0; i != exons.isize(); ++i) {
    int len = min(exons[i].length + pars.JunkOnEnd(),
		  target[exons[i].chromosome].isize() - exons[i].start);
    if (PickBaits(baits, pars, target,
		  exons[i].chromosome, exons[i].start, len,i) > 0)
      ++count;
  }
  return count;
}

