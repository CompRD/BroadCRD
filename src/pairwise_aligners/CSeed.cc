/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "pairwise_aligners/CSeed.h"

/**
 * CSeed
 * Constructor
 */
CSeed::CSeed( ) :
  target_id_ ( -1 ),
  start_ ( 0 ),
  rc_ ( false )
{ }

/**
 * CSeed
 * Constructor
 */
CSeed::CSeed( int target_id, int start, Bool rc ) :
  target_id_ ( target_id ),
  start_ ( start ),
  rc_ ( rc )
{ }

/**
 * CSeed
 * Set
 */
void CSeed::Set( int target_id, int start, Bool rc)
{
  target_id_ = target_id;
  start_ = start;
  rc_ = rc;
}

/**
 * CSeed
 * Print
 */
void CSeed::Print( ostream &out ) const
{
  out << "[t_" << target_id_
      << " @" << start_
      << ( rc_ ? " rc]" : " fw]" );
}

/**
 * CSeed
 * operator<
 */
bool operator< ( const CSeed &left, const CSeed &right ) {
  if ( left.target_id_ == right.target_id_ )
    if ( left.rc_ == right.rc_ )
      return ( left.start_ < right.start_ );
    else return ( left.rc_ < right.rc_ );
  else return ( left.target_id_ < right.target_id_ );
}
