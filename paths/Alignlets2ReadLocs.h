///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__ALIGNLETS_2_READ_LOCS__H
#define PATHS__ALIGNLETS_2_READ_LOCS__H

#include "PairsManager.h"
#include "Vec.h"
#include "VecUtilities.h"
#include "paths/Alignlet.h"
#include "paths/ReadLoc.h"

/**
 * struct SAlignletsData
 *
 * A simple container for core alignment data.
 */
struct SAlignletsData {

  SAlignletsData( const PairsManager *pairs = 0,
		  const vec<alignlet> *aligns = 0,
		  const vec<int> *index = 0 ) :
    pairs_( pairs ), aligns_ ( aligns ), index_ ( index ) { }
  
  void Set( const PairsManager *pairs,
	    const vec<alignlet> *aligns,
	    const vec<int> *index ) {
    pairs_ = pairs; aligns_ = aligns; index_ = index; }

  bool IsEmpty( ) const { return ( pairs_ == 0 ) ; }

  const PairsManager *pairs_;
  const vec<alignlet> *aligns_;
  const vec<int> *index_;
  
};

/**
 * AppendReadLocs
 *
 * Add alignments from given library to the set of read_locs (do
 * nothing if aligns_data is empty).
 *
 * WARNING: locs will not be sorted! This function is the core
 * workhorse for Alignlets2ReadLocs below (which sort locs after all
 * libraries have been added).
 */
void AppendReadLocs( const SAlignletsData &aligns_data,
		     const int readClass,
		     vec<read_loc> &locs,
		     ostream *log = 0 );

/**
 * Alignlets2ReadLocs
 *
 * Convert pairs, alignlets, and matching index file into a (sorted!)
 * vector of read_locs. Not all of frags, jumps, Jumps need to be
 * given (you can pass an empty SAlignletsData for frags, for
 * example).
 */
void Alignlets2ReadLocs( const SAlignletsData &frags,
			 const SAlignletsData &jumps,
			 const SAlignletsData &Jumps,
			 vec<read_loc> &locs,
			 ostream *log = 0 );

#endif
