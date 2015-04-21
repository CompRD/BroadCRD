/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef C_POS_H
#define C_POS_H

#include "system/Types.h"

/**
 * class CPos
 *
 * Base position on a genomic interval. It stores a numeric id, a rc
 * flag, and the position on the interval.
 */
class CPos {

public:

  CPos( ) : id_( -1 ), pos_( -1 ), rc_( False ) { }

  CPos( int id, int pos, Bool rc ) : id_( id ), pos_( pos ), rc_( rc ) { }

  int Id( ) const { return id_; }
  
  int Pos( ) const { return pos_; }
  
  bool IsRc( ) const { return rc_; }

  friend bool operator< ( const CPos &left, const CPos &right ) {
    if ( left.id_ == right.id_ ) {
      if ( left.pos_ == right.pos_ ) {
	return ( left.rc_ < right.rc_ ); }
      return ( left.pos_ < right.pos_ ); }
    return ( left.id_ < right.id_ );
  }
  
  
private:

  int id_;   // numeric id of sequence (for example, read id)
  int pos_;  // position on sequence (for example, position on read)
  Bool rc_;  // tag to specify if pos is on rc of sequence

};

#endif
