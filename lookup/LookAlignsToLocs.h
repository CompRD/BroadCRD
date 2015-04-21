/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef LOOK_ALINGS_TO_LOCS_H
#define LOOK_ALINGS_TO_LOCS_H

#include "ReadLocation.h"
#include "lookup/LookAlign.h"

/**
 * LookAlignsToLocs
 *
 * Generate proto-contigs out of a given set of qlt's.
 *
 * hits: will be sorted
 * locs: the output
 * max_cov: if given, discard contigs with excessive estimated coverage
 * min_reads: if given, discard contigs with less than min_reads reads
 * min_gap: if given, allow some gaps between reads (i.e. bald spots)
 */
void LookAlignsToLocs( vec<look_align_plus> &hits,
		       vec<read_location> &locs,
		       double *max_cov = 0,
		       int *min_reads = 0,
		       int *max_gap = 0,
		       ostream *log = 0 );

/**
 * LookAlignsToGapLocs
 *
 * Generate proto-contigs out of a given set of qlt's. The difference
 * with LookAlignsToLocs above is that here gaps on contigs are
 * allowed. The output consists of one contig per target sequence
 * (exactly). Targets with no hits will be still saved (as contigs
 * with no reads).
 * 
 * tlens: lengths of target sequences
 * hits: will be sorted
 * locs: the output
 */
void LookAlignsToGapLocs( const vec<int> &tlens,
			  vec<look_align_plus> &hits,
			  vec<read_location> &locs,
			  ostream *log = 0 );

#endif
