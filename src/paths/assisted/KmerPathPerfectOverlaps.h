///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__ASSISTED__KMER_PATH_PERFECT_OVERLAPS__H
#define PATHS__ASSISTED__KMER_PATH_PERFECT_OVERLAPS__H

#include "paths/KmerPath.h"
#include "paths/KmerPathInterval.h"

/**
 * GlueKmerPaths
 *
 * It does not test that the two reads overlap, it just glues them at
 * the specified segments.
 */
KmerPath GlueKmerPaths( const KmerPath &left,
			const KmerPath &right,
			const int lpos,
			const int rpos );

/**
 * PerfectOverlap
 *
 * Decide if the kmer_path rpath (the "read") aligns perfectly
 * kmer_path upath (the "unipath"), where the alignment is seeded at
 * segments rpos (on rpath), and upos (on upath).
 */
bool PerfectOverlap( const KmerPath &rpath,
		     const KmerPath &upath,
		     const int rpos,
		     const int upos );

/**
 * IsReadEmbedded
 *
 * Decide if the kmer_path rpath (the "read") is fully embedded in (or
 * exactly equal to) the kmer_path upath (the "unipath"), by a perfect
 * alignment seeded at segments rpos / upos.
 */
bool IsReadEmbedded( const KmerPath &rpath,
		     const KmerPath &upath,
		     const int rpos,
		     const int upos );

/**
 * IsProperExtension
 *
 * Decide if kmer_path rpath (the "read") extends properly kmer_path
 * upath (the "unipath"). This means that the two kmer_paths overlap
 * perfectly by at least one kmer (where the alignment is seeded at
 * segments rpos / upos), and that the read extends properly the
 * unipath on the right.
 */
bool IsProperExtension( const KmerPath &rpath,
			const KmerPath &upath,
			const int rpos,
			const int upos );

/**
 * SimplePrintAlign
 *
 * Simple minded printing utility: unipaths are "glued" at the given
 * pos, one unipath per line. The two unipaths do not need to actually
 * align.
 */
void SimplePrintAlign( const KmerPath &path1,
		       const KmerPath &path2,
		       const int pos1,
		       const int pos2,
		       ostream &out );

#endif
