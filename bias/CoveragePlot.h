/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// CoveragePlot.  Plot coverage from alignments of reads to reference.
//
// Interval to be plotted is TIG.START-STOP.
//
// Coverage = read start points per WINDOW.
//
// GRANULARITY = 1 = default; raise to get finer resolution.
//
// OUT = .png file to write plot to.

#ifndef COVERAGE_PLOT_H
#define COVERAGE_PLOT_H
          
#include "Basevector.h"
#include "CoreTools.h"
#include "lookup/LookAlign.h"

void CoveragePlot( const vecbasevector& reads, const vec<look_align>& aligns,
     const vecbasevector& ref, const int TIG, int START, int STOP,
     const int WINDOW, const int GRANULARITY, const double COVERAGE_DIVIDER,
     const String& OUT );

#endif
