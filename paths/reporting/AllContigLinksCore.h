///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ALL_CONTIG_LINKS_CORE_H
#define ALL_CONTIG_LINKS_CORE_H

#include "Fastavector.h"
#include "PairsManager.h"
#include "Superb.h"
#include "Vec.h"
#include "paths/Alignlet.h"

/**
 * AllContigLinksCore
 *
 * Print all bundles between all contigs.
 *
 * MIN_LINKS: only accept bundles with >= these many links
 * out: output stream
 * log: log stream
 */
void AllContigLinksCore( const int MIN_LINKS,
			 const vec<fastavector> &contigs,
			 const vec<alignlet> &aligns,
			 const vec<int> &index,
			 const PairsManager &pairs,
			 ostream &out,
			 ostream &log );

#endif
