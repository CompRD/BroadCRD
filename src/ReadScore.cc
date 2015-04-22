/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "math/Arith.h"
#include "Qualvector.h"
#include "ReadScore.h"

/**
 * read_score
 * Constructor
 */
read_score::read_score( const base_prob *scorer ) :
  scorer_ ( scorer )
{ }

/**
 * read_score
 * Score
 *
 * For each i call p_i the probability that the base at i is wrong. Then
 * score := 1 - \prod ( 1 - p_i ). In words, this is the probability
 * that at least one base is incorrect.
 */
Float read_score::Score( const qualvector &qual, int start, int len ) const
{
  int end = start + len;
  Float score = 1;
  for (int ii=start; ii<end; ii++)
    score = score * ( Float( 1 ) - scorer_->Prob( qual[ii] ) );
  
  return Float( 1 ) - score;
}

/**
 * read_score
 * Score
 *
 * Same, but quals are now a vec<int>. Since it can be used with tiling
 * type quality scores (in which gaps and missing bases are assigned
 * negative quality scores), it only takes into account the nonnegative
 * entries of the input vector.
 */
Float read_score::Score( const vec<int> &qual, int start, int len ) const
{
  int end = start + len;
  Float score = 1;
  for (int ii=start; ii<end; ii++)
    if ( qual[ii] >= 0 )
      score = score * ( Float( 1 ) - scorer_->Prob( qual[ii] ) );
  
  return Float( 1 ) - score;
}

/**
 * read_score
 * Append
 */
Float read_score::Append( int qual, Float score ) const
{
  ForceAssert( qual >= 0 );
  Float prob = scorer_->Prob( qual );
  return prob + ( score * ( Float( 1 ) - prob ) );
}

