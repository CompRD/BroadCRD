///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Alignment.h"
#include "btl/SanitizeSmithWatBandedA2.h"

/**
 * SanitizeSmithWatBandedA2
 */
bool SanitizeSmithWatBandedA2( align &al )
{
  // Check last block.
  if ( al.Lengths( al.Nblocks( ) -1 ) < 1 ) {
    align loc_al;
    
    // Remove last block, if its length is zero.
    int p1 = al.pos1( );
    int p2 = al.pos2( );
    avector<int> gaps = al.Gaps( );
    avector<int> lens = al.Lengths( );
    int loc_nb = al.Nblocks( ) - 1;
    gaps.resize( loc_nb );
    lens.resize( loc_nb );
    
    // Adjust what is now last length.
    lens( lens.length - 1 ) += Abs( al.Gaps( al.Nblocks( ) - 1 ) );
    
    // Set loc_al sna swap.
    loc_al.Set( p1, p2, gaps, lens );
    swap( loc_al, al );
  }
  
  // Just fail, if there are other bugs in SimithWatBandedA2.
  bool is_valid = true;
  for (int ii=0; ii<al.Nblocks( ); ii++) {
    if ( al.Lengths( ii ) < 1 ) {
      is_valid = false;
      break;
    }
  }

  // Done, return.
  return is_valid;

}

