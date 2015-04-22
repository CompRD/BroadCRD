///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__FIND_VECTOR__H
#define PATHS__FIND_VECTOR__H

#include "Basevector.h"
#include "lookup/LookAlign.h"
#include "pairwise_aligners/KmerAligner.h"

/**
 * FindVectorCore
 *
 * Core aligning function: align one vector to one read. See
 * FindVector below for documentation.
 *
 * v_id: id of vector
 * r_id: id of read
 * v_rc: if vector is rc
 * aligner8: KmerAligner, to find candidate placements
 */
void FindVectorCore( const int v_id,
		     const int r_id,
		     const bool v_rc,
		     const bvec &vbases,
		     const bvec &read,
		     const int MAX_MISMATCHES,
		     const int MAX_INDELS,
		     const double MAX_ER,
		     KmerAligner<8> &aligner8,
		     vec<look_align> &hits );

/**
 * FindVector
 *
 * It looks for vector sequence in a set of reads. This is done by
 * aligning vectors (both fw and rc, or fw only) to the reads, after
 * seeding on 8-mer perfect mathches.
 *
 * HEURISTICS
 *
 * Accept an align if
 *  1) there are <= MAX_MISMATCHES mismatches, and <= MAX_INDELS indels, or
 *  2) the error rate is <= MAX_ER.
 *
 * hits_file: output file (look_aligns)
 * vectors: sequences to screen for
 * reads: set of reads to be screened
 * log: report progress (if not null)
 * hits: if given, store aligns in hits
 */
void FindVector( const String &hits_file,
		 const vecbvec &vectors,
		 const vecbvec &reads,
		 const int MAX_MISMATCHES,
		 const int MAX_INDELS,
		 const double MAX_ER,
		 const bool SKIP_RC,
		 ostream *log = 0,
		 vec<look_align> *hits = 0 );

#endif

