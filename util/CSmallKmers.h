///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef UTIL__C_SMALL_KMERS__H
#define UTIL__C_SMALL_KMERS__H

#include "Basevector.h"
#include "Intvector.h"

/**
 * class CSmallKmers
 *
 * Convert k-mers to int and vice versa (k<=12). Optionally, generate
 * maps to find all k-mers in a given vecbasevector (in which case the
 * two maps target_ids_ and target_pos_ are built - in sync with each
 * other).
 */
class CSmallKmers {
  
public:

  CSmallKmers( const int k, const vecbvec *target = 0 );
  
  // Generate maps for given target.
  void SetTargetBases( const vecbvec &target );
  
  // Value of k.
  int K( ) const { return k_; }

  // Convert k-mer starting at pos in the given bvec.
  int ToInt( const bvec &bases, const int pos ) const;
  
  // Convert -mer of specified target at given pos to int.
  int ToInt( const int tid, const int pos ) const;
  
  // Convert int to k-mer.
  void ToBases( int id, bvec &bases ) const;
  
  // All placements for the int representing a given 8-mer.
  size_t  NPlacements ( const int id )               const;
  int64_t TargetId    ( const int id, const int ii ) const;
  int     TargetPos   ( const int id, const int ii ) const;
  
  
private:

  void Setup( );
  
  
private:
  
  int k_;                   // the k
  vec<int> fourPow_;        // powers of four
  vecbvec target_;          // optional target dataset
  Int64VecVec target_ids_;  // optional map k-mer id to set of target ids
  VecIntVec target_pos_;    // optional map k-mer id to pos on target

};

#endif
