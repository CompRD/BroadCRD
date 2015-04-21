/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef ALIGNS_ON_CONSENSUS_H
#define ALIGNS_ON_CONSENSUS_H

#include "Alignment.h"
#include "Basevector.h"
#include "LocsHandler.h"
#include "Qualvector.h"
#include "tiled/CRealign.h"
#include "tiled/TAlign.h"

/**
 * AlignsOnConsensus
 *
 * It realigns reads to the contig they belong to. If trust_locs is set
 * to true, then it will try to align reads to the place suggested by
 * the given locs file (therefore only aligns respecting the given loc
 * induced orientation will be accepted). Otherwise sort aligns by
 * their score (pick lowest score).
 *
 * Logging levels (beware: serious performance hit if > 1!)
 *  1: print one line per contig,
 *  2: as 1, plus one line per read tested
 *  3: as 2, plus print alignments iff score >= 100.0
 *  4: as 3, plus print all aligmnments
 *
 * the_taligns: output
 * contig_id: the contig id
 * first_locs: position in locs of first read with given contig_id
 * locs: vector of read_location's
 * trust_locs: only accept aligns if they conform with the placements from loc
 * max_discrep: betwen given loc and found align (used iff trust_locs is true)
 * log_level: from 1 (minimum logging) to 4 (maximum logging), see above
 * fast: set low MAXCLIQ and MAX_ALIGNS (faster, but may miss some aligns)
 * max_gap: if not-null, discard an align if it has excessive gaps
 * out: out (log) stream
 */
void AlignsOnConsensus( vec<t_align> &the_taligns,
			int contig_id,
			const vecbasevector &consensus_bases,
			const vecqualvector &consensus_quals,
			const vecbasevector &read_bases,
			const vecqualvector &read_quals,
			const vec<int> &first_locs,
			const vec<read_location> &locs,
			bool trust_locs = false,
			int max_discrep = 128,
			int log_level = 1,
			bool fast = false,
			int *max_gap = 0,
			ostream *out = 0,
			bool checkPerformance = false,
			int numRandomReads = 200 );

/**
 * AlignsOnConsensus
 *
 * A second implementation of AligsOnConsensus: it uses the class
 * CRealign to place reads on the consensus (see the documentation in
 * there for details).
 * 
 * REMARK: setting max_discrep = -1 will accept any placement of the
 * read on the consensus (no check done).
 */
void AlignsOnConsensus( vec<t_align> &the_taligns,
			int contig_id,
			const vecbasevector &cbases,
			const vecbasevector &rbases,
			const lhandler &locs,
			int max_discrep = 48,
			float max_error = 0.35,
			ostream *log = 0,
			ostream *fastalog = 0 );

#endif
 
