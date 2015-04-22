///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef REPORTING__SCAFFOLDED_FASTA_TO_NSTATS_H
#define REPORTING__SCAFFOLDED_FASTA_TO_NSTATS_H

#include "CoverageAnalyzer.h"
#include "Fastavector.h"
#include "math/NStatsTools.h"

/**
 * ScaffoldedFastaToNStats
 *
 * Load sequence from a fasta file, and treat each fasta entry as a
 * scaffolded unit: "N"s (or "n"s) represent either ambiguities or
 * gaps (a gap is a set of >= <GAPN> consecutive "N"s). Report NStats
 * of input scaffolds, contigs, and gaps.
 *
 * GAPN: a gap is a sequence of >= GAPN consecutive "N"s (or "n"s)
 */
void ScaffoldedFastaToNStats( ostream &out,
			      const int GAPN,
			      const vec<fastavector> &fasta );

#endif
