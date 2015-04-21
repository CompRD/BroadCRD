/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef PROGRESSION_H
#define PROGRESSION_H

#include "ParseSet.h"
#include "STLExtensions.h"
#include "Vec.h"

/**
 * GetIntProgression
 *
 * By example. Say seed={2,3,6,10}, min=1000, max=100000. Then,
 * GetIntProgression will return a progression vector {1000, 2000,
 * 3000, 6000, 10000, 20000, 30000, 60000, 100000}.
 *
 * seed: must be sorted, and seed[0] > 1 (if not passed, default will be used).
 */
inline vec<int> GetIntProgression( int min, int max, vec<int> *pseed = 0 )
{
  vec<int> result;

  vec<int> default_seed;
  ParseIntSet( "{2,3,6,10}", default_seed );
  
  const vec<int> &seed = pseed ? *pseed : default_seed;
  if ( pseed ) {
    ForceAssert( is_sorted( seed.begin( ), seed.end( ) ) );
    ForceAssert( seed.size( ) > 0 );
    ForceAssert( seed[0] > 1 );
  }

  int post = min;
  int current = post;
  result.push_back( min );
  while ( current <= max ) {
    for (int ii=0; ii<(int)seed.size( ); ii++) {
      current = seed[ii] * post;
      if ( current > max )
	break;
      result.push_back( current );
    }
    post = current;
  }

  return result;
}

#endif
