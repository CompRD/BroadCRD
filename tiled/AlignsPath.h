// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef ALIGNS_PATH_H
#define ALIGNS_PATH_H

#include "Alignment.h"
#include "ReadLocation.h"
#include "Vec.h"
#include "lookup/LookAlign.h"
#include "tiled/BaggedReads.h"



/*
 * class aligns_path
 *
 * Checks if it is possible to walk from one read to another
 * by jumping only on reads which align well to each other.
 * It is possible to pass null pointers for bases and quals.
 */
class aligns_path {

public:
  
  aligns_path( );

  aligns_path ( const vec<alignment_plus> *aligns,
		const vec<int> *aligns_index,
		const vec<look_align_plus> *hits,
		const bagged_reads *read_bags,
		const vecbasevector *bases = 0,
		const vecqualvector *quals = 0,
		ostream *log = 0 );
  
  void SetPointers( const vec<alignment_plus> *aligns,
		    const vec<int> *aligns_index,
		    const vec<look_align_plus> *hits,
		    const bagged_reads *read_bags,
		    const vecbasevector *bases = 0,
		    const vecqualvector *quals = 0,
		    ostream *log = 0 );
  
  void SetMaxAlignScore( float max_align_score );
  
  void SetMaxPolyScore( int max_poly_score );
  
  void SetUsePolyScore( bool use_poly_score );

  float MaxAlignScore( ) const { return max_align_score_; }
  
  int MaxPolyScore( ) const { return max_poly_score_; }
  
  bool UsePolyScore( ) const { return use_poly_score_; }
  
  bool AreConnected( int pos1, int pos2,
		     const vec<int> &loc_indices,
		     vec<int> &dead ) const;

  bool IsBadLoc( int pos1, const vec<int> &loc_indices ) const;
  
  bool IsDeadForward( int pos1, const vec<int> &loc_indices ) const;
  
  bool IsDeadBackward( int pos1, const vec<int> &loc_indices ) const;

  int CountAligners( int pos1,
		     const vec<int> &loc_indices,
		     vec<int> *aligner_ids = 0 ) const;
  
  
private:

  void Setup( );

  bool LocOverlap( int loc1, int loc2 ) const;

  bool AlignOverlap( int loc1, int loc2 ) const;

  int FindAlign( int loc1, int loc2 ) const;
  
  
private:
  
  const vec<alignment_plus> *aligns_;   // read-read alignments
  const vec<int> *aligns_index_;        // index of the above
  const vec<look_align_plus> *hits_;    // alignments to master sequence
  const bagged_reads *read_bags_;       // original locations
  const vecbasevector *bases_;          // reads bases
  const vecqualvector *quals_;          // reads quals
  ostream *log_;                        // log file (may be null)

  float max_align_score_;               // max value for align score
  int max_poly_score_;                  // max vale for poly score
  bool use_poly_score_;                 // use poly score to salvage aligns
  
};



#endif

