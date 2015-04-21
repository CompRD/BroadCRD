///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PHYSICAL_COVERAGE_FROM_ALIGNS_H
#define PHYSICAL_COVERAGE_FROM_ALIGNS_H

#include "PairsManager.h"
#include "Vec.h"
#include "VecUtilities.h"
#include "paths/TranslateAligns.h"

/**
 * PhysicalCoverageFromAlignsCore
 *
 * Overall physical coverage is computed as the weighted average of
 * the physical coverage of each super (the weight being the length of
 * the super). The physical coverage of a super is defined as N / D,
 * where N and D are defined as follows.
 *
 *  1. An insert is defined "valid" if its two end reads point at each
 *     other, if its stretch does not exceed MAX_STRETCH, and if it is
 *     not closer than RADIUS to the end of the super.
 *
 *  2. For a valid insert, report the length of the portion of the
 *     insert which lies within ( 2 * RADIUS ) bases from the ends of
 *     the super.
 *
 *  3. N is defined as the sum of the reported lengths for all valid
 *     inserts of a super.
 *
 *  4. D is the gapped super length, minus ( 4 * RADIUS ).
 *
 * aligns: these are changed (from contig to super coordinates)
 * MAX_STRETCH: used to define valid inserts
 * RADIUS: see above
 * VERBOSE: print also coverage for each super
 */
double PhysicalCoverageFromAlignsCore( ostream &log,
				       vec<alignlet> &aligns,
				       const vec<int> &index,
				       const vec<superb> &supers,
				       const PairsManager &pairs,
				       const double MAX_STRETCH = 5.0,
				       const int RADIUS = 40000,
				       const bool VERBOSE = true );

#endif
