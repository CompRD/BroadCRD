///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__REPORTING__C_HITS_H
#define PATHS__REPORTING__C_HITS_H

#include "Basevector.h"
#include "CoverageAnalyzer.h"
#include "PairsManager.h"
#include "SeqInterval.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"

/**
 * class CHits
 *
 * Common paired-reads routines: find valid pairs, estimate lib stats,
 * etc. Note that both separations and insert lengths are stored and
 * used in the code (insert length is the separation plus the read
 * lengths).
 */
class CHits {

public:

  // The constructor also estimates lib stats (and fill nds_).
  CHits( const vec<look_align> &hits,
	 const PairsManager &pairs,
	 const vec<int> &cglens,
	 ostream *log = 0 );

  // Find all regions with (physical) coverage for given lib_id <= max_cov.
  void CoveragesAtMost( const int lib_id,
			const int max_cov,
			vec<seq_interval> &regions ) const;
  
  // Find all regions with (physical) coverage for given lib_id = cov.
  void CoveragesExactly( const int lib_id,
			 const int cov,
			 vec<seq_interval> &regions ) const;
  
  // Find all regions with (physical) coverage for given lib_id >= min_cov.
  void CoveragesAtLeast( const int lib_id,
			 const int min_cov,
			 vec<seq_interval> &regions ) const;

  // Update pairs file.
  void UpdatePairsFile( const String &pairs_file, ostream *log = 0 ) const;
  
  
private:
  
  // Clean all data.
  void Cleanup( );
  
  // As usual, to_hit_: (-1 = unmapped), (-2 = multiply placed).
  void GenerateMaps( );
  
  // Fill seps_ (ilens_ is filled only if max_stretch is given).
  void FindSeparations( const double *max_stretch = 0 );

  // Use guess if unable to estimate - HEURISTICS inside, and rather crude!
  void EstimateLibStats( );

  // If align is fully embedded.
  bool IsFullyEmbedded( const look_align &hit ) const;
  
  
private:

  const vec<look_align> &hits_;    // all aligns
  const PairsManager &pairs_;      // pairing info
  const vec<int> &cglens_;         // contig lengths
  ostream *log_;                   // optional log stream
  
  vec<int> to_hit_;                // map read_id to hit_id
  
  vec<NormalDistribution> nds_;   // computed lib stats
  vec< vec<seq_interval> > seps_;  // separations (by library)
  vec< vec<seq_interval> > ilens_; // insert lengths  (by library)

};

#endif
