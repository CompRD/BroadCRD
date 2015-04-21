/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef BASE_PROB_H
#define BASE_PROB_H

#include "math/Arith.h"
#include "Vec.h"

/**
 * \ingroup grp_quals
 *
 * Translates a quality score into the probability that the given base
 * is wrong (if q is the quality score, then the probability is defined as
 * 10^(-q/10)). An inverse operation is provided for convenience.
 *
 * cap_: if true then cap probabilities within a given reasonable bracket.
 *       Min and max values for probabilities are defined in Setup( ).
 *
 * \sa BaseErrorProbProfile::BaseErrorProbProfile(const qualvector&)
 */
class base_prob {
  
public:
  
  base_prob( bool cap = true ) : cap_ ( cap ) { this->Setup( ); }
  
  void SetCap( bool cap ) { cap_ = cap; this->Setup( ); }
  
  inline Float Prob( int qual ) const { return prob_[qual]; }
  
  inline Float Qual( Float prob ) const { return -10 * log10( prob ); }
  
  
private:
  
  void Setup( );
  
  
private:
  
  bool cap_;
  vec<Float> prob_;
  
};

#endif
