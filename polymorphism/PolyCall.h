#ifndef POLY_CALL_H
#define POLY_CALL_H

/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include <iosfwd>

/// \file PolyCall.h
/// A simple struct to hold polymorphism (eg. SNP) info

struct PolyCall {
  int refID;
  int posOnRef;
  unsigned char refBase;  
  unsigned char newBase;  
  PolyCall(int refID_=-1, int posOnRef_=-1, unsigned char refBase_=4,
	   unsigned char newBase_=4) :
    refID(refID_), posOnRef(posOnRef_), refBase(refBase_), newBase(newBase_)
  { }
  /// as_base: convert integer to base (inverse to as_char)
  char as_base(unsigned char x) const;
  /// as_char: convert base to integer (inverse to as_base)
  unsigned char as_char(char x) const;
};

std::ostream & operator<<(std::ostream &out, const PolyCall &p);
std::istream & operator>>(std::istream &in, PolyCall &p);

bool operator<(const PolyCall &a, const PolyCall &b);

#endif //POLY_CALL_H

