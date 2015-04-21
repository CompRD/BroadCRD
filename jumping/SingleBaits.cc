/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "jumping/SingleBaits.h"
#include "math/Functions.h"

void SingleBaitParams::Validate() const {
}

void PickBaits(vec<Bait> & baits, const SingleBaitParams & pars,
	       vecbasevector & target, int contig, 
	       int targetStart, int targetLength, int exonId ) {
  pars.Validate();

  ForceAssertGt(targetLength, pars.JunkOnEnd());

  int baitstart = targetStart;

  while (baitstart < targetStart + targetLength - pars.JunkOnEnd()) {
    int end = baitstart + pars.fragLength;
    Bait site(contig, baitstart, end, exonId);
    baits.push_back(site);
    baitstart = end;
  }
}

void PickAllBaits(vec<Bait> & baits, const SingleBaitParams & pars,
		  vecbasevector & target, const vec<ExonStats> & exons) {
  baits.clear();
  for (int i=0; i != exons.isize(); ++i) {
    int len = min(exons[i].length + pars.JunkOnEnd(),
		  target[exons[i].chromosome].isize() - exons[i].start);
    PickBaits(baits, pars, target, exons[i].chromosome, exons[i].start, len,i);
  }
}
    
