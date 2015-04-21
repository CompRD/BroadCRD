/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef C_TAG_QUAL_H
#define C_TAG_QUAL_H

#include "LocsHandler.h"
#include "Qualvector.h"
#include "SeqInterval.h"

/**
 * class CTagQual
 *
 * Tag regions of contig qualb's (based on several criteria), and save them
 * as seq_intervals. Tag intervals closer than fudge_ are merged together.
 * If tag_edges_ is false, windows at the edges of contigs (begin or end)
 * will not be tgged. Selection criteria may be applied sequentially.
 *
 * Criterion 1: tag all bases with quality < min_good_.
 *
 * Criterion 2: tag all bases for which coverage is < pit_min_cov_. Coverage
 *  is computed by looking at the locs info only, and it takes into account
 *  the "good portions" of reads only. These are defined as follows: walk
 *  along the read bases until there are at least pit_min_size_ bases of
 *  quality >= pit_min_qual_.
 */
class CTagQual {


public:

  CTagQual( const vecqvec *cg_quals = 0 );
  
  // Set quality pointer (contigs).
  void SetCgQuals( const vecqvec *cg_quals );

  // Set quality and locs pointers (reads).
  void SetRPointers( const vecqvec *r_quals, const lhandler *locs );

  void SetFudge( int fudge ) { fudge_ = fudge; }

  void SetTagEdges( bool tag ) { tag_edges_ = tag; }

  // Set min qual for criterion 1.
  void SetMinGood( int min ) { min_good_ = min; }

  // Set "pit" arguments for criterion 2.
  void SetPit( int min_cov, int min_qual, int min_size );

  // Remove tags that do not match any of the given annotations.
  void SelectOverlapping( const vec<seq_interval> & annot );

  // Tag bases with qual < min_good_.
  const vec<seq_interval> &LowQual( bool append = false );
  
  // Tag bases at low coverage (only use good chunks of reads to compute cov).
  const vec<seq_interval> &LowPitCoverage( bool append = false );
  
  // Fetch tags_.
  const vec<seq_interval> &Tags( ) const { return tags_; }

  // Fetch ftags_.
  const vec<int> &Ftags( ) const { return ftags_; }

  // Send info to stream.
  void PrintInfo( ostream &out ) const;
  
  
private:
  
  void SetDefaults( );
  
  void CleanUp( );

  // Generate ftags_ (will sort ftags_).
  void GenerateFtags( );

  // Merge intervals closer than fudge_ (will GenerateFtags at the end).
  void FudgeMerge( );
  
  // Sort, remove windows at edges, etc. In other words, finalize for saving.
  void SortAndRenum( );
  
  
private:

  const vecqvec *cg_quals_; // quals of contigs
  const vecqvec *r_quals_;  // quals of reads (may be null)
  const lhandler *locs_;    // locs data (may be null)

  vec<seq_interval> tags_;  // the result (tag windows, sorted)
  vec<int> ftags_;          // id in tags_ of first tag for given contig

  int fudge_;       // merge low quality windows on contig closer than fudge
  bool tag_edges_;  // do we report low quality regions at ends of contig?

  int min_good_;    // tag regions with qual < min_good_
  
  int pit_min_cov_;  // tag regions at cov < pit_min_cov, where cov is computed
  int pit_min_qual_; //  with good portions of reads (walk on read until at
  int pit_min_size_; //  least pit_min_size bases have qual >= pit_min_qual_

};

#endif
