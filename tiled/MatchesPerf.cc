/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Vec.h"
#include "tiled/MatchesPerf.h"

/**
 * MatchesPerf
 */
bool MatchesPerf( const vec<char> &bases1,
		  const vec<char> &bases2,
		  ostream *plog )
{
  // Log stream.
  ofstream devnull( "/dev/null" );
  ostream &log = plog ? *plog : devnull;

  // Find long and short.
  const vec<char> *b_short = &bases1;
  const vec<char> *b_long = &bases2;
  if ( b_short->size( ) < b_long->size( ) )
    return false;
  if ( b_short->size( ) > b_long->size( ) )
    swap( b_short, b_long );
  
  log << "short_read: |";
  for (int ii=0; ii<(int)b_short->size( ); ii++)
    log << (*b_short)[ii];
  log << "|\nlong_read:  |";
  for (int ii=0; ii<(int)b_long->size( ); ii++)
    log << (*b_long)[ii];
  log << "|\n";

  // Short read has shrunk to nothing.
  if ( b_short->size( ) < 1 ) {
    log << "do not match (short read has 0 length)\n";
    return false;
  }

  // Check identity by attaching the short read on the left.
  bool they_match = true;
  for (int ii=0; ii<(int)b_short->size( ); ii++) {
    if ( (*b_short)[ii] != (*b_long)[ii] ) {
      they_match = false;
      break;
    }
  }
  if ( they_match ) {
    log << "match\n";
    return true;
  }

  // No need to further check.
  if ( b_short->size( ) == b_long->size( ) ) {
    log << "do not match\n";
    return false;
  }
  
  // Check by attaching short on the right.
  they_match = true;
  int short_last = b_short->size( ) - 1;
  int long_last = b_long->size( ) - 1;
  for (int ii=0; ii<(int)b_short->size( ); ii++) {
    if ( (*b_short)[short_last - ii] != (*b_long)[long_last - ii] ) {
      they_match = false;
      break;
    }
  }
  
  log << ( they_match ? "match\n" : "do not match\n" );
  return they_match;
}

