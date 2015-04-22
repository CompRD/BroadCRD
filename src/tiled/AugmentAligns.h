// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef AUGMENT_ALIGNS_H
#define AUGMENT_ALIGNS_H

#include "PackAlign.h"
#include "Basevector.h"
#include "ReadLocation.h"
#include "Qualvector.h"
#include "Vec.h"
#include "VecAlignmentPlus.h"
#include "lookup/LookAlign.h"



/*
 * class augment_aligns
 *
 * The input is basically given in the form of a full set of read-read
 * alignments, and a set of alignments of reads onto a chromosome of
 * another genomic sequence (the master sequence). The idea is to use
 * the information from the latter to decide if two reads are close (and
 * therefore should align) to boost the initial set of read-read aligns.
 * The output is a new set of alignments: the original set (but only
 * for reads aligning the given chromosome), plus the new alignments
 * found as explained above.
 *
 * Options:
 *  copy_only_: do not add any new alignment (just keep the old ones).
 *  max_new_: max number of new alignments tested per read.
 *  max_discrepancy_: max difference between implied and observed overlap.
 *  max_error_rate_: max error rate for an alignment.
 */
class augment_aligns {
  
public:
  
  augment_aligns( );
  
  augment_aligns( const vec<look_align_plus> *hits,
		  const vec_alignment_plus *all_aligns,
		  const vecbasevector *bases,
		  const vecqualvector *quals,
		  const vec<int> *indiv = 0,
		  ostream *log = 0,
		  ostream *new_aligns_log = 0 );
  
  void SetPointers( const vec<look_align_plus> *hits,
		    const vec_alignment_plus *all_aligns,
		    const vecbasevector *bases,
		    const vecqualvector *quals,
		    const  vec<int> *indiv = 0,
		    ostream *log = 0,
		    ostream *new_aligns_log = 0 );
  
  void SetCopyOnly( bool copy ) { copy_only_ = copy; }

  void SetMaxNew( int max ) { max_new_ = max; }
  
  void SetMaxDiscrepancy( int max ) { max_discrepancy_ = max; }
  
  void SetMaxErrorRate( float rate ) { max_error_rate_ = rate; }

  bool CopyOnly( ) const { return copy_only_; }

  int MaxNew( ) const { return max_new_; }

  int MaxDiscrepancy( ) const { return max_discrepancy_; }

  float MaxErrorRate( ) const { return max_error_rate_; }

  void Generate( );
  
  void Generate( int begin_hit, int end_hit );
  
  void SortAndSave( const String &aligns_file );
  
  
private:  
  
  void Setup( );

  void CopyAligns( int begin_hit, int end_hit );
  
  void GenerateNewAligns( int begin_hit, int end_hit );
  
  bool AddAlign( int hit1, int hit2 );

  bool IsInAllAligns( int hit1, int hit2 ) const;
  
  int Overlap( int hit1, int hit2 ) const;

  float AlignBandedSW( int hit1, int hit2, align &al ) const;
  
  float AlignATB( int hit1, int hit2, align &al ) const;
  
  float ErrorRate( int hit1, int hit2, const align &al ) const;
  
  
private:

  const vec<look_align_plus> *hits_;
  const vec_alignment_plus *all_aligns_;
  const vecbasevector *bases_;
  const vecqualvector *quals_;
  const vec<int> *indiv_;    // vector with individual ids (may be null)
  ostream *log_;             // basic log (may be null)
  ostream *new_aligns_log_;  // new alignments log (may be null)
  
  bool copy_only_;           // do not add any new alignment
  int max_new_;              // max number of reads tested for align per read
  int max_discrepancy_;      // between implied and observed overlap amount
  float max_error_rate_;     // max error rate in a read-read alignment
  mutable ofstream devnull_;

  vec<alignment_plus> sub_aligns_;

};



#endif
