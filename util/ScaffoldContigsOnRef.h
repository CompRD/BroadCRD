/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef UTIL__SCAFFOLD_CONTIGS_ON_REF__H
#define UTIL__SCAFFOLD_CONTIGS_ON_REF__H

#include "Basevector.h"
#include "Fastavector.h"
#include "Superb.h"
#include "lookup/LookAlign.h"
#include "paths/SaveScaffoldGraph.h"

/**
 * ScaffoldContigsOnRef
 *
 * Use alignments of contigs onto a reference to generate a scaffold
 * structure (no links needed). initial_supers are used to tag which
 * contigs need to be stored in the output scaffolds: if a contig
 * appears in initial_supers but it is not aligned to the reference,
 * then it will appear in output as a singleton super.
 *
 * MAX_GAP: max gap allowed between adjacent contigs
 * MAX_OVERLAP: max overlap allowed between adjacent contigs
 * MIN_GAP_DEV: min threshold for gap dev (usually set to gap_size / 4 )
 */
void ScaffoldContigsOnRef( const vec<superb> &initial_supers,
			   const vec<fastavector> &initial_contigs,
			   const vec<look_align> &aligns,
			   const String &ASSEMBLY_OUT,
			   const int MAX_GAP,
			   const int MAX_OVERLAP,
			   const int MIN_GAP_DEV,
			   ostream *log = 0 );

#endif

