///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__ASSISTED__C_CLOSURES__H
#define PATHS__ASSISTED__C_CLOSURES__H

#include "Basevector.h"

/**
 * enum closure_type
 *
 * Description codes for closures (or partial closures) from a cloud
 * walk.
 */
enum closure_t {
  OVERLAP,        // overlap
  FULL_CLOSURE,   // full closure
  LEFT_CLOSURE,   // partial closure (left side)
  RIGHT_CLOSURE   // partial closure (right side)
};

/**
 * class CClosures
 *
 * Container for the data structure needed to represent all closures
 * between two contigs.  The three core member vectors have the same
 * size (and are in sync).  Each closure can be one (and only one) of
 * the following
 *
 * An overlap:
 *   bases_[ii].size( ) = 0,
 *   overlaps_[ii] >= 1,
 *   ctypes_[ii] = OVERLAP.
 *
 * A true closure:
 *   bases_[ii].size( ) >= 0,
 *   overlaps_[ii] = -1,
 *   ctypes_[ii] != OVERLAP.
 *
 */
class CClosures {

public:

  CClosures( ) { this->Clear( ); }
  
  // Clear, reserve, etc.
  void Clear( );
  void Reserve( size_t nclosures );
  void Resize( size_t nclosures );

  // Core adding methods.
  void AddOverlap( const int over );
  void AddTrueClosure( const bvec &bases, const closure_t ct );
  void AddClosure( const bvec &bases, const int over, const closure_t ct );
  void AddClosure( const CClosures &other, const int ii );
  
  // Remove (K-1) bases from left and/or right of closures.
  void ChompAdjacencies( const int K );

  // Full closures and overlaps first (sorted by distance from gap_size).
  void GapBasedSort( int gap_size );
  
  // Const accessors.
  size_t Size( ) const { return bases_.size( ); }
  bool IsOverlap( int ii ) const { return ctypes_[ii] == OVERLAP; }
  bool IsFull( int ii ) const { return ctypes_[ii] == FULL_CLOSURE; }
  bool IsLeft ( int ii ) const { return ctypes_[ii] == LEFT_CLOSURE; }
  bool IsRight ( int ii ) const { return ctypes_[ii] == RIGHT_CLOSURE; }
  int Overlap( int ii ) const { return overlaps_[ii]; }
  int Len( int ii ) const { return bases_[ii].size( ); }
  closure_t Ctype( int ii ) const { return ctypes_[ii]; }
  const bvec &Bases( int ii ) const { return bases_[ii]; }
  const vecbvec &Bases( ) const { return bases_; }
  const vec<int> &Overlaps( ) const { return overlaps_; }
  const vec<closure_t> &Ctypes( ) const { return ctypes_; }
  
  // Info line (for specified event).
  void PrintInfo( int ii, int super_id, int gap_id, ostream &out ) const;

  // Simple report (one liner).
  void PrintReport( ostream &out ) const;
  
  
private:
  
  vecbvec bases_;
  vec<int> overlaps_;
  vec<closure_t> ctypes_;
  
};

#endif
