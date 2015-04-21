/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef ALIGN_PAIRS_TO_HYPER_H
#define ALIGN_PAIRS_TO_HYPER_H

#include "Basevector.h"
#include "Bitvector.h"
#include "CoreTools.h"
#include "PairsManager.h"
#include "feudal/BinaryStreamTraits.h"
#include "graph/Digraph.h"
#include "paths/HyperBasevector.h"
#include "paths/HyperKmerPath.h"
#include "paths/SeqOnHyper.h"


// Computes alignments of reads to the hyperkmer path directly using Perfectlookup

typedef pair<int,int> IndexPair;

void AlignPairsToHyper( const HyperKmerPath& h, const HyperBasevector& hb,
			// note that hb must be HyperBasevector(h)
			const vecbasevector& reads, const PairsManager& pairs,
			const String& sub_dir,
			vec<CompressedSeqOnHyper>& csaligns,
			vec<vec<int> >& csaligns_index,
			vec< IndexPair >& pair_places, 
			Bool verbose,  const int max_insert_size = -1,
                        vecbasevector* p_cached_edges = 0,
                        vec<partial_perfect_align>* p_cached_ppaligns = 0);

// EvaluatePair: given alignments ap1, ap2 of the reads of a read pair P to a
// HyperKmerPath, determine if the insert could be placed contiguously on the
// HyperKmerPath, without stretching too much.  Return True in that case.

Bool EvaluatePair( 
     const SeqOnHyper& ap1,            // alignment of first read
     const SeqOnHyper& ap2,            // alignment of second read
     const int& sep, const int& stdev, // the read pair
     double dev_mult,                  // how much stretch is allowed in devs
     const HyperBasevector& hb,        // from the HyperKmerPath
     const vec<int>& to_right_vertex,  // index from edge to right vertex
     const vec<int>& to_left_vertex,   // index from edge to left vertex
     const digraphE<int>& G,           // mirrors hb but edges are lengths in kmers
     const vec<IndexPair>& paired,
     const vec<Bool>& have_path,
     const vec<int>& min_path,
     const vec<int>& max_path,
     const vec<int>& a,
     const vec<int>& b );

// EvaluatePair: given alignments cs1, cs2 of the reads of a read pair P to a
// HyperKmerPath, determine if the insert could be placed contiguously on the
// on a single edge of the HyperKmerPath, without stretching too much.
// Return True in that case.
// Special case of the more general EvalutePair function. Optimized for single
// edge condition using CompressedSeqOnHyper information only.

Bool EvaluatePairOnSingleEdge( 
     const CompressedSeqOnHyper& cs1,     // alignment of first read
     const CompressedSeqOnHyper& cs2,     // alignment of second read
     const int& sep, const int& stdev,    // the read pair
     double dev_mult );                   // how much stretch is allowed in devs

// Find placements of read pairs that are not too stretched.  The returned data is
// in pair_places, and we guarantee that csaligns[ pair_places[i].first ] is
// forward aligned, whereas csaligns[ pair_places[i].second ] is reverse aligned.
// Note that pairs of alignments going between components in the HyperKmerPath are 
// not included.

void FindPairPlacements( const HyperKmerPath& h, const HyperBasevector& hb,
			 const vec<CompressedSeqOnHyper>& csaligns,
			 const vec<vec<int> >& csaligns_index,
			 const int nreads, const PairsManager& pairs,
			 vec< IndexPair >& pair_places, Bool verbose,
			 const String& sub_dir, const int max_insert_size = -1,
			 const int timeout = 100, Bool skip_on_timeout = False);


// Find the bases on the HyperBasevector edges which are covered by alignments of 
// read pairs that are not too stretched.  Note that read pairs going between 
// components are not counted (which may be a good thing or a bad thing).  One can
// filter to only use pairs within a given separation range and to only use 
// forward- or reverse-aligned reads.

void FindWellCovered( int min_sep, int max_sep, Bool use_fw, Bool use_rc,
     const HyperKmerPath& h, const HyperBasevector& hb,
     const vec<CompressedSeqOnHyper>& csaligns, 
     const vec< IndexPair >& pair_places,
     int nreads, const PairsManager& pairs,
     vecbitvector& cov, Bool verbose );


// Deletes edges of HyperKmerPath that are not covered by bases in a good pair.

void FindPoorlyCovered( HyperKmerPath& h, HyperBasevector& hb,
     const KmerBaseBroker& kbb, const vec<CompressedSeqOnHyper>& csaligns, 
     const vec< IndexPair >& pair_places, const vecbasevector& reads, 
     const PairsManager& pairs, Bool remove, int verbosity,
     const Bool NEW_POORLY_COVERED = True );

/**
   Class: Bridge

   Represents a link between two edges of one contig of the assembly,
   comprised of a set of read pairs with one read of the pair landing on one edge
   and the other read on the other edge.  This class stores the fact that there
   is a bridge, and summary information about the set of pairs comprising the bridge:
   the separation they imply between the edges (mean and stddev), number of pairs,
   number of pairs where one read is unique, number of pairs where both reads are
   unique.

   See also: <WritePairAligns>
*/
class Bridge {
 public:
  int edge1, edge2;
  nbases_t sepMean, sepStd;
  int npairs, npairsOneReadUnique, npairsBothReadsUnique;

  Bridge() { }
  Bridge( int _edge1, int _edge2, nbases_t _sepMean, nbases_t _sepStd,
	  int _npairs, int _npairsOneReadUnique, int _npairsBothReadsUnique ):
    edge1(_edge1), edge2(_edge2), sepMean(_sepMean), sepStd(_sepStd),
    npairs(_npairs), npairsOneReadUnique(_npairsOneReadUnique), npairsBothReadsUnique(_npairsBothReadsUnique) {
  }

};  // class Bridge
TRIVIALLY_SERIALIZABLE(Bridge);

#endif
