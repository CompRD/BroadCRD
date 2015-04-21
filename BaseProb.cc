/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "math/Arith.h"
#include "BaseProb.h"
#include "math/Functions.h"

/**
 * base_prob
 * Setup
 *
 * Build a lookup table for the probability of calling a wrong base.
 */
void base_prob::Setup( )
{
  int prob_size = 256;
  Float min_prob = 0.75;      // used only if cap_ = true
  Float max_prob = 0.00001;   // used only if cap_ = true
    
  prob_.resize( prob_size );

  for (int ii=0; ii<prob_size; ii++) {
    prob_[ii] = Pow( Float( 10 ), - ( Float( ii ) / Float( 10 ) ) );
    if ( cap_ ) {
      prob_[ii] = Min( min_prob, prob_[ii] );
      prob_[ii] = Max( max_prob, prob_[ii] );
    }
  }  
}
