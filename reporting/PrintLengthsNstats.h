/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef PRINT_LENGTHS_NSTATS_H
#define PRINT_LENGTHS_NSTATS_H

#include "PrettyPrintTable.h"
#include "String.h"
#include "Vec.h"
#include "math/Functions.h"

/**
 * PrintLengthsNstats
 *
 * Print a simple text-based N-distribution plot for the given vector
 * of lengths. It will generate output in a form similar to that of
 * BasicAssemblyStats.
 */
void PrintLengthsNstats( const vec<int> &lens,
			 const String &name,
			 ostream &out );

#endif
