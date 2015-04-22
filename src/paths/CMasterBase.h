///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__C_MASTER_BASE__H
#define PATHS__C_MASTER_BASE__H

#include "Basevector.h"

/**
 * class CMasterBase
 *
 * It encapsulates all the variants of a base on a given sequence.
 * For each event for this base (see below), there are maps to the ids
 * of objects (reads or pairs) supporting that event. Events can be a
 * letter base (A, C, G, T), a deletion, or an insertion of one ore
 * more bases.
 *
 * HEURISTICS
 *
 * min_covs_: interval of min coverages for an haplotype to be
 *   considered (it is given as an interval because we require a
 *   larger min coverage for bases at higher coverage).
 */
class CMasterBase {

public:
  
  CMasterBase( );
  
  // Clear all data structures.
  void ClearAll( );

  // Heuristics are documented in SetDefauls( ), in the .cc file.
  void SetMinCovs( pair<double,double> min_covs );

  // Add event.
  void AddBase      ( const size_t &id, const char base );
  void AddInsertion ( const size_t &id, const bvec &ins );
  void AddDeletion  ( const size_t &id );

  // Counters.
  size_t InsCoverage ( ) const;
  size_t DelCoverage ( ) const;
  size_t BaseCoverage( ) const;
  size_t TotCoverage ( ) const;

  // Most common base as a char (use as_base( ) to convert to A, C, G, or T).
  char MostCommonBase( ) const;

  // Print a brief descriptor (one line).
  void PrintBrief( const String &info, ostream &out ) const;
  
  // Valid haplotypes (ie events with cov>= min_coverage_).
  void ListHaplotypes(const String &info, const char base, ostream &log) const;
  int NHaplotypes( ) const;
  
  
private:
  
  void SetDefaults( );

  // Compute required min coverage (based on min_covs_ and total coverage).
  double MinCoverage( ) const;
  
  
private:

  pair<double,double> min_covs_;   // heuristics - see Defaults( ) in the .cc
  mutable double min_coverage_;    // see ComputeCoverage( ) in the .cc
  
  vec<size_t> idsA_;
  vec<size_t> idsC_;
  vec<size_t> idsG_;
  vec<size_t> idsT_;
  vec<size_t> idsDels_;
  vec< vec<size_t> > idsIns_;
  vecbvec basesIns_;
  
};

#endif
