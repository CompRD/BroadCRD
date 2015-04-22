/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Vec.h"
#include "tiled/CharBaseUtil.h"
#include "tiled/PadsRatio.h"

/**
 * PadsRatio
 */
float PadsRatio( const vec<char> &consensus, const vec<char> &read )
{
  ForceAssert( read.size( ) == consensus.size( ) );
  ForceAssert( read.size( ) > 0 );

  int consensus_len = 0;
  int read_pads = 0;
  for (int ii=0; ii<(int)consensus.size( ); ii++) {
    if ( IsGap( consensus[ii] ) ) continue;
    consensus_len++;
    if ( IsGap( read[ii] ) ) read_pads++;
  }
  if ( consensus_len < 1 ) return 1.0;
  
  return SafeQuotient( read_pads, consensus_len );
}

