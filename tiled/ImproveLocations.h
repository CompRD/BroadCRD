// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef IMPROVE_LOCATIONS_H
#define IMPROVE_LOCATIONS_H

#include "Alignment.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "Vec.h"
#include "VecAlignmentPlus.h"
#include "tiled/BaggedReads.h"



/*
 * class improve_locations
 *
 * Given an initial set of locations for a protocontig, it tries to
 * fine tune up the StartOnContig's by using read-read alignments,
 * and to single out (and eliminate) reads which do not seem to
 * align consistently with adjacent reads.
 */
class improve_locations {
  
public:
  
  improve_locations( );
  
  improve_locations( const vec_alignment_plus *aligns,
		     const vec<read_location> *locs,
		     ostream *log = 0,
		     ostream *exclude_log = 0 );
  
  void SetPointers( const vec_alignment_plus *aligns,
		    const vec<read_location> *locs,
		    ostream *log = 0,
		    ostream *exclude_log = 0 );
  
  void Improve( );
  
  void Improve( int contig_id );
  
  void SortAndFill( vec<read_location> &improved_locs );
  
  void SortAndSave( const String &out_file );
  
  
private:
  
  void Setup( );
  
  void EliminateBadLocs( );

  void AdjustLocs( );

  void OverlapByLoc( int id, vec<int> &mini_ids ) const;

  void OverlapLeftByLoc( int id, vec<int> &mini_ids ) const;
  
  void OverlapRightByLoc( int id, vec<int> &mini_ids ) const;
  
  int GetAlignIndex( int mini1, int mini2 ) const;

  int LocOverlap( int mini1, int mini2 ) const;

  int LocOffset( int min1, int min2 ) const;

  int AlignOverlap( int mini1, int mini2 ) const;

  int AlignOffset( int mini1, int mini2 ) const;
  
  bool EvalShift( int mini1, int mini2, int &shift ) const;

  
private:
  
  const vec_alignment_plus *aligns_;  // read-read aligns
  const vec<read_location> *locs_;    // original locations
  ostream *log_;                      // log file (may be null)
  ostream *exclude_log_;              // excluded reads log (may be null)
  
  vec<read_location> ilocs_;          // improved locations
  vec<int> first_locs_;               // first locs for original locations
  ofstream devnull_;

  vec<read_location> minilocs_;       // temporary locs for building contigs
  
};



#endif
