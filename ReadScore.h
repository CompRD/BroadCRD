/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef READ_SCORE_H
#define READ_SCORE_H

#include "math/Arith.h"
#include "BaseProb.h"
#include "Qualvector.h"

/**
 * class read_score
 *
 * It assignes a score to any stretch of a given read, based only on its
 * quality scores. The higher the scores of a fixed stretch, the lower
 * its score will be. Range: the closed interval [0, 1].
 */
class read_score {
  
public:
  
  read_score( const base_prob *scorer );
  
  Float Score( const qualvector &qual, int start, int len ) const;
  
  Float Score( const vec<int> &qual, int start, int len ) const;

  // Append qual to given score and return new score.
  Float Append( int qual, Float score ) const;
  
  
private:
  
  const base_prob *scorer_;
  
};

#endif
