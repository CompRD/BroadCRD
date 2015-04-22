///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Alignment.h"
#include "btl/WalkAlign.h"

/**
 * WalkOn1
 */
pair<int,int> WalkOn1( const align &al, const int on1 )
{
  const int nan = numeric_limits<int>::max( );

  // Out of range.
  if ( on1 < al.pos1( ) ) return make_pair( -1, 0 );
  if ( on1 > al.Pos1( ) ) return make_pair( nan, 0 );

  // Special case, on1 = pos1.
  if ( on1 == al.pos1( ) ) return make_pair( al.pos2( ), 0 );

  // Walk on align.
  int p1 = al.pos1( );
  int p2 = al.pos2( );
  for (int block=0; block<al.Nblocks( ); block++) {
    int gap = al.Gaps( block );
    if ( gap > 0 ) 
      p2 += gap;
    if ( gap < 0 ) {
      for (int ii=0; ii<-gap; ii++) {
	p1++;
	if ( p1 == on1 ) return make_pair( p2, -gap - ii );
      }
    }
    int len = al.Lengths( block );
    for (int ii=0; ii<len; ii++) {
      p1++;
      p2++;
      if ( p1 == on1 ) return make_pair( p2, 0 );
    }
  }
  
  // We should not get here.
  ForceAssert( 1 == 0 );
  return make_pair( 666, 666);
  
}

/**
 * WalkOn2
 */
pair<int,int> WalkOn2( const align &al, const int on2 )
{
  const int nan = numeric_limits<int>::max( );

  // Out of range.
  if ( on2 < al.pos2( ) ) return make_pair( -1, 0 );
  if ( on2 > al.Pos2( ) ) return make_pair( nan, 0 );

  // Special case, on2 = pos2.
  if ( on2 == al.pos2( ) ) return make_pair( al.pos1( ), 0 );

  // Walk on align.
  int p1 = al.pos1( );
  int p2 = al.pos2( );
  for (int block=0; block<al.Nblocks( ); block++) {
    int gap = al.Gaps( block );
    if ( gap > 0 ) {
      for (int ii=0; ii<gap; ii++) {
	p2++;
	if ( p2 == on2 ) return make_pair( p1, gap - ii );
      }
    }
    if ( gap < 0 )
      p1 += -gap;
    int len = al.Lengths( block );
    for (int ii=0; ii<len; ii++) {
      p1++;
      p2++;
      if ( p2 == on2 ) return make_pair( p1, 0 );
    }
  }
  
  // We should not get here.
  ForceAssert( 1 == 0 );
  return make_pair( 666, 666);
  
}
