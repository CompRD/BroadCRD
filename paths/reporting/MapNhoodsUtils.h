///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__REPORTING__MAP_NHOODS_UTILS_H
#define PATHS__REPORTING__MAP_NHOODS_UTILS_H

// MakeDepend: dependency MakeLookupTable
#include "lookup/LookAlign.h"
#include "paths/reporting/ReftigUtils.h"
#include "util/RunCommand.h"

/**
 * EstimatePlacement
 *
 * Guess placement of local assembly based on the alignment of its
 * seed. It returns: (target_id, begin, end).
 */
triple<int,int,int> EstimatePlacement( const HyperBasevector &hyper,
				       const look_align &la,
				       const vec<int> &alens );

/**
 * RefinePlacement
 *
 * Align the given local assembly on the segment of assembly
 * identified by the given placement. If no estimated placement is
 * given, then align the local assembly against the whole
 * assembly. Return (-1, -1, -1 ) if a nhood cannot be mapped.
 *
 * placement: both input and output
 * aligns: alignments of local nhood unibases onto assembly
 */
void RefinePlacement( triple<int,int,int> &placement,
		      vec<look_align> &aligns,
		      const int K,
		      const int nhood_id,
		      const String &assembly_head,
		      const String &outdir_head,
		      const vec<fastavector> &assembly_fasta,
		      const vec<HyperBasevector> &hypers,
		      const Bool FORCE );

/**
 * MapRange
 *
 * Parse given RANGE and map local nhoods to it. RANGE must be of the
 * form SUPER_ID:BEGIN-END.
 */
void MapRange( ostream &log,
	       const String RANGE,
	       const Bool SKIP_UNMAPPED,
	       const Bool FORCE,
	       const int K,
	       const int n_seeds,
	       const String &out_dir,
	       const String &assembly_head,
	       const vec<fastavector> &assembly_fasta,
	       const vec<int> &assembly_lengths, 
	       const vec<int> &seeds_ids,
	       const vec<int> &aligns_index,
	       const vec<look_align> &seeds_aligns,
	       const vec<HyperBasevector> &hypers );

#endif
