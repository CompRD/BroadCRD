// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
// 

#include "Alignment.h"
#include "Basevector.h"
#include "JustifyPads.h"



/*
 * JustifyPadsRight
 */
void JustifyPadsRight( align &al,
		       const basevector &bases1,
		       const basevector &bases2 )
{
  // pos1: cursor on read1, pos2: cursor on read2.
  int block_id = 0;
  int pos1 = al.pos1( );
  int pos2 = al.pos2( );
  
  // There are no gaps.
  if ( al.Nblocks( ) < 2 )
    return;
  
  // Loop over align's blocks.
  while ( block_id < al.Nblocks( ) - 1 ) {
    int length = al.Lengths( block_id );
    
    // Move to next gap.
    pos1 += length;
    pos2 += length;
    
    // Shift pads and move on.
    while ( pos1 < static_cast<int>(bases1.size( )) &&
	    pos2 < static_cast<int>(bases2.size( )) &&
	    bases1[ pos1 ] == bases2[ pos2 ] ) {
      al.AddToLength( block_id, 1 );
      al.AddToLength( block_id + 1, -1 );
      al.Compactify( bases1.size( ), bases2.size( ) );
      pos1++;
      pos2++;
    }
    int gap = al.Gaps( block_id + 1 );
    if ( gap > 0 )
      pos2 += gap;
    else
      pos1 += -gap;

    // Next block.
    block_id++;
  }
}



