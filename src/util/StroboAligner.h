///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef STROBO_ALIGNERRE_H
#define STROBO_ALIGNERRE_H

#include "Basevector.h"
#include "STLExtensions.h"
#include "util/SearchFastb2Core.h"

/**
 * StroboAligner
 *
 * Cheap semi-perfect aligner.  Break the query sequences into non
 * overlapping chunks, and align perfectly the first K bases of each
 * chunk.  If enough chunks align consistently, then an alignment is
 * declared by the implicit placements of the chunks (warning! This
 * will not be tested with a real aligner!)
 *
 * Remark: at this time only fw alignments are found.
 *
 * aligns: the output (in the same format as SearchFastb2)
 * TARGET: fastb of target sequence
 * QUERY: fastb of query sequence
 * OUTDIR: directory for output files
 * K: select the first K bases from each interval
 * STROBE: size of consecutive intervals
 * MAX_GAP: used to cluster consistent placements
 * MIN_LENGTH: do not even try to align short queries
 * MIN_RATIO: how many consistent alignments are needed to align a query
 */
void StroboAligner( vec< triple<int64_t,int64_t,int> > &aligns,
		    const String &TARGET,
		    const String &QUERY,
		    const String &OUTDIR,
		    const int K = 40,
		    const int STROBE = 100,
		    const int MAX_GAP = 24,
		    const int MIN_LENGTH = 300,
		    const double MIN_RATIO = 0.95 );

#endif
