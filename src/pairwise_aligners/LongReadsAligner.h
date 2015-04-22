///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PAIRWISE_ALIGNERS__LONG_READS_ALIGNER__H
#define PAIRWISE_ALIGNERS__LONG_READS_ALIGNER__H

#include "Basevector.h"
#include "PrintAlignment.h"
#include "lookup/LookAlign.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "util/CSmallKmers.h"

/**
 * SRange (range for a placement of a query onto a target). NB: weight
 * is signed (sign captures orientation).
 */
struct SRange {

  SRange( int tid, int start, int stop, int weight ) :
    target_ ( tid ), start_ ( start ), stop_ ( stop ), weight_ ( weight ) { }

  void PrintInfo( const int query_id, ostream &out ) const {
    out << " q." << query_id
	<< ( weight_ < 0 ? "[-]" : "[+]" )
	<< "  vs  t." << target_ << " at [" << start_ << ", " << stop_ << "]"
	<< "  w = " << Abs( weight_ )
	<< "\n";
  }
  
  friend bool operator< ( const SRange &left, const SRange &right ) {
    if ( Abs( left.weight_ ) > Abs( right.weight_ ) ) return true;
    if ( Abs( left.weight_ ) < Abs( right.weight_ ) ) return false;
    if ( left.target_ < right.target_ ) return true;
    if ( left.target_ > right.target_ ) return false;
    if ( left.start_ < right.start_ ) return true;
    if ( left.start_ > right.start_ ) return false;
    return ( left.stop_ < right.stop_ );
  }
  
  int target_;
  int start_;
  int stop_;
  int weight_;

};



/**
 * LongReadsCoverages
 *
 * Generate a table of coverages from the given aligns. Note:
 *
 * 1. It extrapolates the coverage based on sample_size and total
 *    target length, if the aligns were generated on a random subset.
 * 2. If a query owns more than one alignment, then only the longest
 *    alignment will be counted.
 *
 * csv: if true, generate csv output (comma separated entries)
 * nickname: used only if csv is true (optional nick name for line)
 */
void LongReadsCoverages(const vecbvec &targets,
                        const vecbvec &queries,
                        const vec<look_align> &aligns,
                        const int sample_size,
                        ostream &out,
                        const bool csv = false,
                        const String nickname = "" );

/**
 * LongReadsAligner
 *
 * It aligns queries onto targets by first seeding via shared kmers
 * (kmer size specified by smallk), and then aligning with
 * Smith-Waterman.
 *
 * qltout_file: output file (aligns sorted by query_id)
 * qselect: optionally only align these query ids (does not need to be sorted)
 * VERBOSE: optional verbose log (sent to cout: from 0=no_log, to 2=max_log)
 * log: optional log stream, to show progress
 */
void LongReadsAligner(const int smallk,
                      const vecbvec &targets,
                      const vecbvec &queries,
                      vec<look_align> * look_aligns_p,
                      const vec<int> *qselect = 0,
                      const int VERBOSE = 0,
                      ostream *log = 0 );

#endif
