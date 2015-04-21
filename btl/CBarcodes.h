///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef BTL__C_BARCODES_H
#define BTL__C_BARCODES_H

#include "Intvector.h"
#include "PairsManager.h"
#include "TokenizeString.h"
#include "btl/Minial.h"

/**
 * enum delpair_t
 *
 * Type of deletion for a pair.
 */
enum delpair_t {
  
  NOT_DELETED,    // not deleted
  CAP_LOW_QUAL,   // one or both reads have too many low qual bases
  CAP_UNPLACED,   // cap not placed
  CAP_MULT,       // cap placed twice on the same read
  CAP_ILLOGICAL,  // cap placed on both reads, with illogical orientation
  IDX_MISMATCH,   // cap on both reads, but inconsistent indexes
  IDX_A_CUT,      // not enough space for index A, before cap
  IDX_A_MISS,     // there is space for A, but A is not found
  IDX_B_CUT,      // not enough space for B, before A
  IDX_B_MISS,     // there is space for B, but B is not found
  IDX_C_CUT,      // not enough space for C, before B
  IDX_C_MISS,     // there is space for C, but C is not found
  IDX_D_CUT,      // not enough space for D, before C
  IDX_D_MISS,     // there is space for D, but D is not found
  TOO_SHORT,      // not enough is left after removing cap
  IS_DUPLICATE    // this pair is a duplicate of another  
  
};



/**
 * class CBarcodes
 *
 * A set of four barcode indexes, with the list of all pair ids
 * sharing that set of indexes.
 *
 */
class CBarcodes {

public:

  CBarcodes( ) { this->Clear( ); }

  // Set methods.

  void Clear( );
  void SetOne( const int ii, const int index );
  void SetIndexes( const IntVec &indexes );
  void SetIndexes( const String &as_string );
  void SetAll( const int64_t pid, const IntVec &indexes );
  void Reserve( const size_t n_pids );
  void AddPairId( const int64_t pid );
  
  // Const methods.

  int Rank( ) const;
  size_t size( ) const { return indexes_.size( ); }
  size_t NPairs( ) const { return pids_.size( ); }
  const int &operator[]( size_t ii ) const { return indexes_[ii]; }
  const int64_t &PairId( size_t ii ) const { return pids_[ii]; }
  String CompactName( ) const;
  void Write( ostream &out ) const;
  void PrintOneLine( ostream &out ) const;

  // Friends.

  friend bool operator== ( const CBarcodes &left, const CBarcodes &right );
  friend bool operator!= ( const CBarcodes &left, const CBarcodes &right );
  friend bool operator< ( const CBarcodes &left, const CBarcodes &right );
  
  
private:

  IntVec indexes_;   // four indexes (A, B, C, and D)
  Int64Vec pids_;    // ids of pairs indexed by indexes_

};



/**
 * IntVecRank
 *
 * The rank is defined as the longest uninterrupted stretch of
 * non-negative indexes. For example: ( a b c d ) has rank 4; ( a b c
 * -1 ) has rank 3; and ( a -1 c d ) has rank 1. All barcodes with -1
 * at the A position have rank 0.
 */
int IntVecRank( const IntVec &indexes );

/**
 * IntVecsConsistent
 *
 * Check if two reads from the same pair are assigned the same index.
 */
bool IntVecsConsistent( const IntVec &left, const IntVec &right );

/**
 * WriteVecBarcodes
 */
void WriteVecBarcodes( const String &outfile, const vec<CBarcodes> &barcodes );

/**
 * LoadVecBarcodes
 */
void LoadVecBarcodes( const String &infile, vec<CBarcodes> &barcodes );

/**
 * BuildVecBarcodes
 *
 * Package raw barcodes into a vec<CBarcodes>. It assumes that
 * raw_barcodes and exp_ranks are in sync with the pairs file, and
 * that raw_barcodes[ii].pids_ = {ii}.
 */
void BuildVecBarcodes( const VecIntVec &raw_barcodes,
		       const IntVec &exp_ranks,
		       vec<CBarcodes> &barcodes );

/**
 * DiscardNonFullRank
 *
 * Tag as deleted pairs for which no full-rank barcode is found. Same
 * assumptions as BuildVecBarcodes above.
 */
void DiscardNonFullRank( const int FULL_RANK,
			 const VecIntVec &raw_barcodes,
			 const IntVec &exp_ranks,
			 vec<int> &deleted );
/**
 * SetCapEnds
 * 
 * CapEnd[read_id] contains the base position after the cap ends at each read_id
 * if cap is not found, or the cap aligns as RC, CapEnd[read_id] takes on zero
 * 
 */
void
SetCapEnds(const vecbvec&reads
          ,const bvec& cap_bases
          ,const vec<minial>& caps
          ,vec<short>& CapEnds
          );

/**
 * DiscardShortGenomic
 *
 * Discard all pairs for which the genomic portion is too small, and
 * store in klens the lengths of the portions of the reads containing
 * the genomic portion, up to the start of the cap:
 *
 *              klen
 *           <------->
 *           ---------[rc_cap][...]->
 *
 * klen may be equal to the read length, if the mate of the read
 * contains the cap, and the read does not contain the cap.
 *
 * REMARK: only full-rank pairs (i.e. pairs with all four indexes
 * placed) are considered in input.
 *
 * REMARK: you must make sure that the cap is placed consistently on
 * at least one read of each non deleted pair (in other words, you
 * must run DiscardNonFullRank before DiscardShortGenomic).
 */
void DiscardShortGenomic( const int MIN_GENLEN,
			  const PairsManager &pairs,
			  const vecbvec &reads,
			  const vec<minial> &caps,
			  const vec<CBarcodes> &barcodes,
			  vec<short> &vCapEnds,
			  vec<int> &klens,
			  vec<int> &deleted );

/**
 * DiscardDuplicates
 *
 * Remove duplicates.
 *
 * REMARK: you need to run DiscardShortGenomic before Dedup.
 */
void DiscardDuplicates( const PairsManager &pairs,
			const vecbvec &reads,
			const vecqvec &quals,
			const vec<CBarcodes> &barcodes,
			const vec<short> &vCapEnds,
			const vec<int> &klens,
			vec<int> &deleted );

/**
 * SaveIndexedfasta
 *
 * Save genomic chunks, encoding in the fasta names indexes and
 * original read ids. It saves fasta and quals (plus some helper
 * files).
 */
void SaveIndexedFasta( const PairsManager &pairs,
		       const vecbvec &reads,
		       const vecqvec &quals,
		       const String &out_head,
		       const vec<CBarcodes> &barcodes,
		       const vec<short> &vCapEnds,
		       const vec<int> &klens,
		       const vec<int> &deleted );

/**
 * ReportDeleted
 */
void ReportDeleted( const vec<int> &deleted,
		    const PairsManager &pairs,
		    ostream &out );

#endif
