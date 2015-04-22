///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ALIGN_SEQS_TO_HYPER_H
#define ALIGN_SEQS_TO_HYPER_H

#include "Alignment.h"
#include "Basevector.h"
#include "Intvector.h"
#include "CoreTools.h"
#include "ReadLocationLG.h"
#include "lookup/LookAlign.h"
#include "paths/HyperBasevector.h"
#include "paths/HyperKmerPath.h"
#include "paths/SeqOnHyper.h"

// Core Alignment Code - intended to be called directly:

// Finds all prefect alignments to the HyperKmerPath - PerfectLookup Method (Fastest)
// A sequence must align entirely along the HyperKmerPath - partial alignments are not
// found. Note that hb must be HyperBasevector(h)
void AlignSeqsToHyper( const HyperKmerPath& h, const HyperBasevector& hb,
		       const vecbasevector& seqs, const String& sub_dir,
		       vec<CompressedSeqOnHyper>& csaligns, Bool verbose,
                       vecbasevector* p_cachedEdges = 0,
                       vec<partial_perfect_align>* p_cachedAligns = 0 );

// Finds alignments to the HyperBasevector - PerfectAligner Method (Slow)
// Computes alignments using the PerfectAligner and extends them along the 
// the HyperBasevector.
// This version returns the aligns in the SeqOnHyper form
void AlignSeqsToHyper( const HyperBasevector& hb, const vecbasevector& seqs,
     vec<SeqOnHyper>& saligns );

// Finds alignments to the HyperBasevector - PerfectAligner Method (Slow)
// Computes alignments using the PerfectAligner and extends them along the 
// the HyperBasevector.
// This version returns the aligns in the CompressedSeqOnHyper form
void AlignSeqsToHyper( const HyperBasevector& hb, const vecbasevector& seqs,
     vec<CompressedSeqOnHyper>& csaligns );

// Finds alignments to the HyperBasevector - Alternative Method (Faster)
// Uses precomputed read alignments to unipaths to determine alignments to
// the HyperBasevector. This should be quicker than the PrefectAligner based
// methods of calculating alignments.
void AlignSeqsToHyper( const HyperKmerPath& h, const HyperBasevector& hb,
     const vecbasevector& seqs, const vecbasevector& unibases, 
     const vecKmerPath& unipaths, const vec<ReadLocationLG>& ulocs, 
     const VecULongVec& ulocs_indexr, const String& sub_dir, const String& tmprun,
     vec<alignment_plus>& aligns );


// Aux Alignment Code - not intended to be called directly:


// Finds alignments to the HyperBasevector - Extend alignment_plus Method
// Extends previously found alignments along the HyperBasevector.
// Requires a Vec of alignment_plus objects that contain alignments to each
// edge of the HyperBasevector. It extends these alignments along the 
// HyperBasevector (where possible), returning the result as either a Vec
// of SeqOnHyper or CompressedSeqOnHyper, or both.
void AlignSeqsToHyper( const HyperBasevector& hb, const vecbasevector& seqs,
     const vec<alignment_plus>& aligns, 
     vec<SeqOnHyper>* p_saligns, vec<CompressedSeqOnHyper>* p_csaligns );


// Wrapper around AlignSeqToHyper (Extend alignment_plus Method)
// This version returns the aligns in the SeqOnHyper from only
// Called by the PerfectAligner based version of AlignSeqToHyper
void AlignSeqsToHyper( const HyperBasevector& hb, const vecbasevector& seqs,
     const vec<alignment_plus>& aligns, vec<SeqOnHyper>& saligns );


// Wrapper around AlignSeqToHyper (Extend alignment_plus Method)
// This version returns the aligns in the CompressedSeqOnHyper from only
// Called by the PerfectAligner based version of AlignSeqToHyper
void AlignSeqsToHyper( const HyperBasevector& hb, const vecbasevector& seqs,
     const vec<alignment_plus>& aligns, vec<CompressedSeqOnHyper>& csaligns );


// Misc Alignment Code:


// Compresses alignments, extending any partial alignments along the hyperkmerpath.
// Called by the PerfectLookup based version of AlignSeqToHyper
void ExtendAndCompressAlignments ( const HyperKmerPath& h, const HyperBasevector& hb,
				   const vecbasevector& reads, const vecbasevector& edges,
				   const vec<partial_perfect_align>& aligns,
				   vec<CompressedSeqOnHyper>& csaligns );

// Builds an index into the alignments based on the alignment Id1 value
void BuildAlignmentIndex (const vec<CompressedSeqOnHyper>& csaligns, int nreads, 
			  vec<vec<int> >& csaligns_index, const bool verbose = false );

// Useful only for debugging purposes.  Returns #erros -- we hope 0.
int VerifyAlignSeqsToHyper( const HyperBasevector& hb, 
			    const vecbasevector& seqs );


#endif
