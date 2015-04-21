///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__MAP_ALIGNLET_TO_REF_H
#define PATHS__MAP_ALIGNLET_TO_REF_H

#include "Vec.h"
#include "lookup/LookAlign.h"
#include "paths/Alignlet.h"

/**
 * MapAlignletToRef
 *
 * Given the alignment of a read on a contig, and of contigs onto a
 * reference, convert the alignment of the read in reference
 * coordinates (multiple alignments are returned, if a contig aligns
 * the reference on more than one place). No alignment is returned if
 * the contig does not align the reference.
 *
 * INPUT
 *  al_onC: alignlet of read on contig
 *  cg_hits: alignments of contigs on reference
 *  to_hits: map contig_id to ids of alignments of contig_id on reference

 * OUTPUT:
 *  als_onR: alignlets on reference (possibly empty)
 *  cg_hit_ids: if not null, specify ids of cg_hits (in sync with als_onR)
 */
void MapAlignletToRef( const alignlet &al_onC,
		       const vec<look_align_plus> &cg_hits,
		       const vec< vec<int> > &to_hits,
		       vec<alignlet> &als_onR,
		       vec<int> *cg_hit_ids = 0 );

#endif
