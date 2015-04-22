///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// AssessAssemblyChangesCore.
// Experimental.

#ifndef ASSESS_ASSEMBLY_CHANGES_CORE_H_
#define ASSESS_ASSEMBLY_CHANGES_CORE_H_

#include "Fastavector.h"
#include "CoreTools.h"

void assess_assembly_changes(const String & TARGETS, 
                             const String & scaffolds_tigs_file,
                             const String & efasta_file,
                             const String & fasta_file,
                             const String & data_dir,
                             const vec<int> & tigs);


#endif
