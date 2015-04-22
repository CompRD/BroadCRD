/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef FOUR_BASE_2_H
#define FOUR_BASE_2_H

#include "Basevector.h"
#include "math/Functions.h"
#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"

class four_base2 {

public:

  four_base2( ) { }
  four_base2( float A, float C, float G, float T )
  {    base_[0] = A;
  base_[1] = C;
  base_[2] = G;
  base_[3] = T;   }

  float A( ) const { return base_[0]; }
  float C( ) const { return base_[1]; }
  float G( ) const { return base_[2]; }
  float T( ) const { return base_[3]; }

  float base(int i) { return base_[i]; }

  friend ostream& operator<<( ostream& out, const four_base2& b );

  // Call: call a base.  This always returns an answer, in the case of equality
  // favoring A > C > G > T.  Return 0, 1, 2, or 3.

  int Call( ) const;

  // CallQuality: return a primitive metric with value between 1 and infinity
  // (1000000000).  It's the ratio of the value for the best base to the next
  // best base.  The worse value (1) should be interpreted for now as 
  // "completely random".

  float CallQuality( ) const;

  float MaxInt() const { return *max_element(base_,base_+4); }

  float Sum() const { return accumulate(base_,base_+4, 0.f); }
  
  friend
  bool operator==(const four_base2 & lhs, const four_base2 & rhs) {
    return 0 == memcmp(lhs.base_, rhs.base_, 4 * sizeof(float));
  }

  bool ShortOk() const { return *max_element(base_, base_+4) < 32768
		     && *min_element(base_, base_+4) > -32768; }

private:

  float base_[4];

};

TRIVIALLY_SERIALIZABLE(four_base2);
typedef SerfVec<four_base2> FourBase2Vec;
typedef MasterVec<FourBase2Vec> VecFourBase2Vec;
extern template class SmallVec< four_base2, MempoolAllocator<four_base2> >;
extern template class OuterVec<FourBase2Vec>;

inline bool operator !=(const four_base2 & lhs, const four_base2 & rhs) {
  return !(lhs == rhs);
}

/// Call bases for an entire intensity vector:
void Call( const VecFourBase2Vec& I, vecbasevector& bases );

#endif //FOUR_BASE_2_H
