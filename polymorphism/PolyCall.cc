/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "polymorphism/PolyCall.h"
#include "CoreTools.h"
#include "Basevector.h"
#include <iostream>

/// as_base: convert integer to base (inverse to as_char)
char PolyCall::as_base(unsigned char x) const
{
  if (x<4)
    return ::as_base(x);
  return '-';
}

/// as_char: convert base to integer (inverse to as_base)
unsigned char PolyCall::as_char(char x) const
{
  if (x=='-')
    return 4;
  return ::as_char(x);
}
  


ostream & operator<<(ostream &out, const PolyCall &p)
{
  out << p.refID << " " << p.posOnRef << " ";
  out << p.as_base(p.refBase) << " " << p.as_base(p.newBase);
  return out;
}

istream & operator>>(istream &in, PolyCall &p)
{
  PolyCall tmp;
  in >> tmp.refID >> tmp.posOnRef;
  string tmpR, tmpN;
  in >> tmpR >> tmpN;
  if (in) {
    ForceAssertEq(tmpR.size(),1u);
    ForceAssertEq(tmpN.size(),1u);
    tmp.refBase = p.as_char(tmpR[0]);
    tmp.newBase = p.as_char(tmpN[0]);
    p=tmp;
  }
  return in;
}
    
bool operator<(const PolyCall &a, const PolyCall &b)
{
  if (a.refID < b.refID) return true;
  if (a.refID > b.refID) return false;
  if (a.posOnRef < b.posOnRef) return true;
  if (a.posOnRef > b.posOnRef) return false;
  if (a.newBase < b.newBase) return true;
  if (a.newBase > b.newBase) return false;
  return false;
}
