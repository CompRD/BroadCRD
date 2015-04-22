///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSEMBLY_INTEGRITY_H
#define ASSEMBLY_INTEGRITY_H

#include "Basevector.h"
#include "Fastavector.h"
#include "Superb.h"
#include "Vec.h"
#include "math/Functions.h"

/**
 * AssemblyIntegrity
 *
 * Run a series of tests to formally validate an allpaths assembly
 * (returns the number of errors found). Only tests for which data are
 * passed are run:
 *
 * if ( fastas && fastbs ): test that fasta and fastb have the same
 *   size, and that the entries sizes are in sync.
 *
 * if ( superbs && ( fastas || fastbs ) ): check contig lengths match
 *   the length given in the superb structure (if both fastbs and
 *   fastbs are given, the contig lengths are determined by the
 *   fastbs).
 *
 * if ( fastbs && efastas ): test that the fastbs match exactly the
 *   efastas (it assumes the efastas are loaded with LoadEfastaFlat)
 */
int AssemblyIntegrity( ostream *log = 0,
		       const vec<superb> *superbs = 0,
		       const vec<fastavector> *fastas = 0,
		       const vecbasevector *fastbs = 0,
		       const vecbasevector *efastas = 0 );

#endif
