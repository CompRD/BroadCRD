///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__ASSISTED__S_PLACEM__H
#define PATHS__ASSISTED__S_PLACEM__H

#include "Basevector.h"
#include "PairsManager.h"
#include "SupersHandler.h"
#include "VecUtilities.h"
#include "paths/ReadLoc.h"

/**
 * struct SPlacem
 *
 * Placement of a read onto a super.
 */
struct SPlacem {

  SPlacem( );

  SPlacem( int64_t rid, int32_t cid, int sid, int beg, int end, bool rc );

  int FindAlign( const vec< triple<int64_t,int64_t,int> > &aligns ) const;
  
  bool IsLogicalPair( const int MIN_SEP,
		      const int MAX_SEP,
		      const SPlacem &other ) const;

  void AddReadLoc( const vec< triple<int64_t,int64_t,int> > &aligns,
		   const vec<int> &clens,
		   const PairsManager &pairs,
		   vec<read_loc> &locs ) const;
  
  void AddReadLocs( const vec< triple<int64_t,int64_t,int> > &aligns,
		    const vec<int> &clens,
		    const PairsManager &pairs,
		    const SPlacem &partner,
		    vec<read_loc> &locs ) const;
  
  friend bool operator< ( const SPlacem &left, const SPlacem &right )
  {
    if ( left.rid_ < right.rid_ ) return true;
    if ( left.rid_ > right.rid_ ) return false;
    // NB: skipping cid_.
    if ( left.sid_ < right.sid_ ) return true;
    if ( left.sid_ > right.sid_ ) return false;
    if ( left.beg_ < right.beg_ ) return true;
    if ( left.beg_ > right.beg_ ) return false;
    if ( left.end_ < right.end_ ) return true;
    if ( left.end_ > right.end_ ) return false;
    return ( left.rc_ < right.rc_ );
  }
  
  int64_t rid_;   // read id
  int32_t cid_;   // contig id (signed)
  int sid_;       // super id
  int beg_;       // begin on super
  int end_;       // end on super
  bool rc_;       // rc on super

};

#endif
