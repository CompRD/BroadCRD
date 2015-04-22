///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "paths/reporting/ReftigUtils.h"
#include <omp.h>
// MakeDepend: library OMP

/**
 * RunGetAlignsFast
 *
 * Call the function GetAlignsFast in ReftigUtils.h (see documentation
 * in there). WARNING: notice that GetAlignFast is omp-parallelized!
 *
 * NUM_THREADS: use all if 0
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_Int( K );
  CommandArgument_String( FASTB_FILE );
  CommandArgument_String( LOOKUP_FILE );
  CommandArgument_String( ALIGNS_FILE );
  CommandArgument_String( TMP_DIR );
  CommandArgument_Bool_OrDefault( USE_CACHE, True );
  CommandArgument_Int_OrDefault( NUM_THREADS, 0 );
  EndCommandArguments;

  NUM_THREADS = configNumThreads( NUM_THREADS );
  omp_set_num_threads( NUM_THREADS );

  vec<look_align> aligns;
  GetAlignsFast( K, FASTB_FILE, LOOKUP_FILE, ALIGNS_FILE, aligns,
		 USE_CACHE, TMP_DIR );

  cout << Date( ) << ": " << aligns.size( ) << " aligns found" << endl;
  
}
