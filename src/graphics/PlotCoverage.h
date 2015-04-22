/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// Utilities to plot hopping window (see below) coverage.

#ifndef PLOT_COVERAGE_H
#define PLOT_COVERAGE_H

#include "Qualvector.h"
#include "SeqInterval.h"
#include "String.h"
#include "Vec.h"

/**
 * PlotCoverage
 *
 * An hopping window is almost like a sliding window, with the difference
 * that instead of sliding one base at a time we slide by ( window_size / 2 )
 * bases at a time (hence the name: we hop from one base to onother with
 * jumps of magnitude ( window_size / 2 ). PlotCoverage plots the hopping
 * window coverage of a genomic interval. Input consists of a vector of
 * floats (the coverages). Output is saved as an eps file.
 */
void PlotCoverage( int window_size,
		   const vec<float> &cov,
		   const String &outfile );

/**
 * PlotCoverage
 *
 * As above, but the input consists of a vector of seq_intervals generated,
 * for instance, by CoverageAnalyzer::GetAllCoverages( ). It returns false
 * if the plot cannot be generated.
 *
 * coverages: must be sorted
 * fcov: first entry in coverages for the given target_id
 */
bool PlotCoverage( int target_id,
		   int target_len,
		   int window_size,
		   const vec<seq_interval> &coverages,
		   const vec<int> &fcov,
		   const String &eps_file );

/**
 * PlotCoverage
 *
 * As above, but the input is given as a qualvector (here coverage is
 * capped at 255). You can use GenerateCoverageVecvec to generate the
 * vecqualvector of coverages.
 */
bool PlotCoverage( int window_size,
		   const qualvector &cov,
		   const String &eps_file );

#endif  // PLOT_COVERAGE_H
