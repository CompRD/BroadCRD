#ifndef BASE_MAP_INFO_H
#define BASE_MAP_INFO_H
/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "Basevector.h"
#include "String.h"
#include "Vec.h"
#include "TokenizeString.h"

struct BaseMapInfo {

  char ref;
  int pos;
  unsigned short rawF[4], rawR[4], accF[4], accR[4], insF[4], insR[4], delF, delR;
  unsigned short allAcc, allRaw;
  unsigned short called[4], inscalled[4], delcalled;

  BaseMapInfo(String & s) {
    vec<String> x;
    Tokenize(s, x);
    Set(x);
  }

  BaseMapInfo(vec<String> &x) { Set(x); }

  BaseMapInfo()
  {
    ref = pos = delF = delR = allAcc = allRaw = delcalled = 0;
    for (int base=0; base<4; ++base) {
      rawF[base] = rawR[base] = accF[base] = accR[base] = insF[base] = insR[base] = called[base] = inscalled[base] = 0;
    }
  }

  void Set(vec<String> &x)
  {
    static const unsigned int NO_INDELS = 20;
    static const unsigned int WITH_INDELS = 32;
    Assert(x.size() == NO_INDELS || x.size() == WITH_INDELS);
    ref = as_char(x[1][0]);
    pos = x[0].Int();

    int base;
    for (base=0; base<4; ++base) {
      // extract bases read from line
      rawF[base] = x[2+base].Int();
      rawR[base] = x[6+base].Int();
      accF[base] = x[10+base].Int();
      accR[base] = x[14+base].Int();
      called[base] = 0;
      inscalled[base] = 0;
    }
    delcalled = 0;

    if (x.size() == WITH_INDELS) {
      for (base=0; base<4; ++base) {
	insF[base] = x[20+base].Int();
	insR[base] = x[24+base].Int();
      }
      delF = x[29].Int();
      delR = x[30].Int();
    }
    else {
      for (base=0; base<4; ++base) {
	insR[base] = insF[base] = 0;
      }
      delR = delF = 0;
    }

    // We omit insertions from sum because they are compatible with
    // other calls, e.g. if 100 reads agree w/ reference and all 100
    // call insertion as well that is 100% of reads agreeing with
    // insertion, not 50%.
    allAcc = delF + delR;
    // We include deletions for now in the raw sum even though that
    // represents only accepted deletion reads.  The map format
    // would have to change to record the raw deletion reads.
    allRaw = delF + delR;

    for (base=0; base<4; ++base) {
      allAcc += accF[base] + accR[base];
      if (accF[base]>0 && accR[base]>0) {
	++called[base];
      }
      if (insF[base]>0 && insR[base]>0) {
	++inscalled[base];
      }
      allRaw += rawF[base] + rawR[base];
    }
    if (delF>0 && delR>0) {
      ++delcalled;
    }
  }

  unsigned short acceptF() const {
    unsigned short sum=0;
    for (int i=0; i != 4; ++i) {
      sum+= accF[i];
    }
    return sum;
  }

  unsigned short acceptR() const {
    unsigned short sum=0;
    for (int i=0; i != 4; ++i) {
      sum+= accR[i];
    }
    return sum;
  }

  unsigned int snipRawF() const {
    return rawF[(ref+1)%4] + rawF[(ref+2)%4] + rawF[(ref+3)%4];
  }

  bool anySnipF() const {
    return accF[(ref+1)%4] || accF[(ref+2)%4] || accF[(ref+3)%4];
  }
  
  unsigned int snipAccF() const {
    return accF[(ref+1)%4] + accF[(ref+2)%4] + accF[(ref+3)%4];
  }

  bool anyInsF() const {
    for (int i=0; i != 4; ++i) {
      if (insF[i]) return true;
    }
    return false;
  }

  bool anyPolyF() const {
    return anySnipF() || anyInsF() || delF;
  }

  char whichSnipF() const {
    for (int i=1; i != 4; ++i) {
      if (accF[(ref+i) % 4]) return i;
    }
    return -1;
  }

  char whichInsF() const {
    for (int i=0; i != 4; ++i) {
      if (insF[i]) return i;
    }
    return -1;
  }

  unsigned int snipRawR() const {
    return rawR[(ref+1)%4] + rawR[(ref+2)%4] + rawR[(ref+3)%4];
  }

  bool anySnipR() const {
    return accR[(ref+1)%4] || accR[(ref+2)%4] || accR[(ref+3)%4];
  }
  
  bool anyInsR() const {
    for (int i=0; i != 4; ++i) {
      if (insR[i]) return true;
    }
    return false;
  }

  bool anyPolyR() const {
    return anySnipR() || anyInsR() || delR;
  }

  char whichSnipR() const {
    for (int i=1; i != 4; ++i) {
      if (accR[(ref+i) % 4]) return i;
    }
    return -1;
  }

  char whichInsR() const {
    for (int i=0; i != 4; ++i) {
      if (insR[i]) return i;
    }
    return -1;
  }

  unsigned int snipAccR() const {
    return accR[(ref+1)%4] + accR[(ref+2)%4] + accR[(ref+3)%4];
  }

  unsigned int snipAcc() const {
    return snipAccR() + snipAccF();
  }

  bool anySnip() const { return anySnipF() || anySnipR(); }
  bool anyIns() const { return anyInsF() || anyInsR(); }
  bool anyDel() const { return delF || delR; }
  bool anyPoly() const { return anyPolyF() || anyPolyR(); }
  
    
  bool Called(int i, bool fwrc = true, int minReads = 1, double minRatio = 0) const {
    unsigned short reads = accR[i] + accF[i];
    return (reads >= minReads)
      && (!fwrc || called[i])
      && (0==minRatio || double(reads)/double(allAcc) >= minRatio);
  }

  bool InsCalled(int i, bool fwrc = true, int minReads = 1, double minRatio = 0) const {
    unsigned short reads = insR[i] + insF[i];
    return (reads >= minReads)
      && (!fwrc || inscalled[i])
      && (0==minRatio || double(reads)/double(allAcc) >= minRatio);
  }

  bool DelCalled(bool fwrc = true, int minReads = 1, double minRatio = 0) const {
    unsigned short reads = delF + delR;
    return (reads >= minReads)
      && (!fwrc || delcalled)
      && (0==minRatio || double(reads)/double(allAcc) >= minRatio);
  }

  double PercentCalled(int i) const {
    return 100.*double(accF[i]+accR[i])/double(allAcc);
  }
  double PercentInsCalled(int i) const {
    return 100.*double(insR[i] + insF[i])/double(allAcc);
  }
  double PercentDelCalled() const {
    return 100.*double(delR + delF)/double(allAcc);
  }

};

bool operator<(const BaseMapInfo & l, const BaseMapInfo & r);

ostream & operator<<(ostream &out, const BaseMapInfo &b);
#endif
