/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef C_SEED_H
#define C_SEED_H

#include "ReadLocation.h"
#include "String.h"
#include "STLExtensions.h"

/**
 * CSeed
 *
 * A seed represents a possible way to align a read to a target.
 */
class CSeed {

public:

  CSeed( );

  CSeed( int target_id, int start, Bool rc );

  void Set( int target_id, int start, Bool rc);
  
  int TargetId( ) const { return target_id_; }
  
  int Start( ) const { return start_; }
  
  Bool RC( ) const { return rc_; }
  
  void Print( ostream &out ) const;
  
  friend bool operator< ( const CSeed &left, const CSeed &right );
  
  
private:

  int target_id_;  // id of target of seed
  int start_;      // start of read on target implied by seed
  Bool rc_;        // orientation of read on target implied by seed

};

#endif
