/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "polymorphism/BaseMapInfo.h"

bool operator<(const BaseMapInfo & l, const BaseMapInfo & r)
{
  return l.pos < r.pos;
}

ostream & operator<<(ostream &out, const BaseMapInfo &b)
{
  const int accF = (b.accF[0]+b.accF[1]+b.accF[2]+b.accF[3]);
  const int accR = (b.accR[0]+b.accR[1]+b.accR[2]+b.accR[3]);
  // Yes, printf really does simplify life, alas. We use sprintf
  // because RunTime() delinks c and c++ stdout's, and the n variant
  // to ensure safety. 
  const int bufsize = 256;
  char buffer[bufsize];
  const char summary_format[] = "%3d/%3d acc (%3d/%3d rej):  ";
  snprintf(buffer, bufsize, summary_format, accF, accR,
	   (b.rawF[0]+b.rawF[1]+b.rawF[2]+b.rawF[3]-accF),
	   (b.rawR[0]+b.rawR[1]+b.rawR[2]+b.rawR[3]-accR));
  out << buffer;
  const char base_format[] = "%c: %3d/%3d (%3d/%3d) ";
  int base;
  for (base=0; base<4; ++base) {
    snprintf(buffer, bufsize, base_format, as_base(base), b.accF[base], b.accR[base], 
	     (b.rawF[base]-b.accF[base]), (b.rawR[base]-b.accR[base]));
    out << buffer;
  }
  const char ins_format[] = "i%c: %3d/%3d ";
  for (base=0; base<4; ++base) {
    snprintf(buffer, bufsize, ins_format, as_base(base), b.insF[base], b.insR[base]);
    out << buffer;
  }
  out << " Del: " << b.delF << "/" << b.delR << " ";
  return out;
}
