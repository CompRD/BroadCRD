///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


/* EvalUtilsLG: Tools for evaluating a HyperKmerPath in the RunAllPathsLG
 * pipeline.  These tools are adapted from the code in EvalHyper.cc and
 * EvalUtils.cc.
 *
 * All functions in this class are titled "Eval...".
 *
 * Josh Burton
 * November 2009
 *
 ******************************************************************************/



#include "Basevector.h"
#include "Bitvector.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "lookup/LookAlign.h"
#include "paths/EvalUtils.h"
#include "paths/HyperKmerPath.h"



// Write a summary of this assembly to output_file.
// Summary includes graph statistics and a listing of gaps.
// If genome and aligns are provided, summary also includes coverage info
// and a number of errors.
// If genome_diploid is provided, the genome is taken as diploid.
void
EvalSummary( const String & output_file,
	     const HyperKmerPath & hkp,
	     const vecbasevector& genome, const vecbasevector& genome_diploid,
	     const vecbitvector& genome_amb,
	     const vec<look_align> & aligns,
	     const vec< vec<int> > & aligns_index );


// Prints library stats info to output_file.
// Wrapper for PairsManager::printLibraryStats.
void
EvalLibraryStats( const String & output_file, const PairsManager & pairs );



// Prints scaffolding accuracy information to output_file.
// Wrapper for the module ScaffoldAccuracy.
void
EvalScaffoldAccuracy( const String & output_file,
		      const String & genome_file,
		      const String & scaffold_file );
