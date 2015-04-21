///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PAIRWISE_ALIGNERS__C_REF_MERGER_H
#define PAIRWISE_ALIGNERS__C_REF_MERGER_H

#include "Basevector.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "pairwise_aligners/SmithWatBandedA.h"

/**
 * class CRefMerger
 *
 * Merge contigs together, based on the overlaps implied by the
 * alignments of the contigs onto a reference. Input aligns and
 * contigs are digested, and contigs are recursively merged together.
 *
 * Remark: the initial set of aligns must be clean (no rc aligns, and
 * only proper aligns).
 */
class CRefMerger {

public:
  
  CRefMerger( const int min_overlap,
	      const int swband_ratio,
	      const int max_gap,
	      const int min_gap,
	      const int min_gap_dev,
	      const vecbvec &ref );

  CRefMerger( const int min_overlap,
	      const int swband_ratio,
	      const int max_gap,
	      const int min_gap,
	      const int min_gap_dev,
	      const vecbvec &ref,
	      const vecbvec &in_contigs,
	      const vec<look_align> &in_aligns );
  
  // Digest original contigs aligns.
  void Setup( const vecbvec &in_contigs, const vec<look_align> &in_aligns );
  
  // Merge contigs: run MergeIterations until no more merges are found.
  void Merge( ostream *log = 0 );

  // Save contigs (only keep contigs >= min_clen).
  void Save( const String &out_head, const int min_clen, ostream *log = 0 );
  
  // Report count and N50 of merged contigs.
  String BriefStats( const int iter, const int *merges = 0 ) const;
  
  
private:

  // One iteration of merging.
  int MergeIteration( );

  // Merge contigs (and aligns of contigs to reference).
  bool MergeContigs( const int idx1, const int idx2 );
  
  // Remove empty contigs/aligns (only keep contigs >= min_clen).
  void CompactifyData( const int min_clen = 0 );
  
  // Find pos on 1 corresponding to pos2 on 2 (-1 on error).
  int PosOn1( const int pos2, const align &al ) const;

  // Find pos on 2 corresponding to pos1 on 1 (-1 on error).
  int PosOn2( const int pos1, const align &al ) const;

  // If align is perfect (no mismatches/indels).
  bool IsPerfect( const bvec &b1, const bvec &b2, const align &al ) const;

  // Overlap implied by aligns on reference.
  int ImpliedOverlap( const int idx1, const int idx2 ) const;
  
  // Internal consistency test for given align.
  bool IsValid( const int idx, ostream *log = 0 ) const;
  
  // Generate a printable string representing the align.
  String BriefAlignInfo( const align &al ) const;
  
  
private:

  const int min_overlap_;    // min overlap to merge contigs
  const int swband_ratio_;   // sw band := ( overlap / swband_ratio )
  const int max_gap_;        // max gap size allowed
  const int min_gap_;        // min gap size allowed
  const int min_gap_dev_;    // min value for gaps' dev

  const vecbvec ref_;        // fastb of target
  vecbvec orcontigs_;        // original contigs
  vecbvec orcontigsU_;       // original contigs (unmapped)
  vec<look_align> oraligns_; // aligns of original contigs

  vecbvec contigs_;          // merged contigs
  vec<look_align> aligns_;   // aligns of merged contigs
  vec< vec<int> > orig_ids_; // ids of original contigs in each merged contig
  
};

#endif
