///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef EXTRACT_FERMI
#define EXTRACT_FERMI

// Extract a fermi assembly.

#include "CoreTools.h"
#include "other_assemblers/ConnectionAssembly.h"

void ExtractSGA( const String& IN_DIR, connection_assembly& A, bool bVariants );

#endif
