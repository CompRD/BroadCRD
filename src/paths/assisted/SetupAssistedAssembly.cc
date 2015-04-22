///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "util/RunCommand.h"
// MakeDepend: dependency RunAllPathsLG

/**
 * SetupAssistedAssembly
 *
 * Prepare data for assisting.
 */
int main( int argc, char *argv[] )
{
  RunTime( );
 
  BeginCommandArguments;

  // AllPathsLG standard args.
  CommandArgument_String( PRE );
  CommandArgument_String( REFERENCE_NAME );
  CommandArgument_String( DATA_SUBDIR );
  CommandArgument_String( RUN );
  CommandArgument_Int_OrDefault( K, 96 );
  CommandArgument_Bool_OrDefault( OVERWRITE, True );
  CommandArgument_UnsignedInt_OrDefault( NUM_THREADS, 0 );

  EndCommandArguments;
  
  // Run AllPathsLG with specified targets.
  String strK = ToString( K );

  String targets_run
    = String( "\"{" )
    + String( "all_reads.unipaths.k" + strK )
    + String( ",all_reads.unipathsdb.k" + strK )
    + String( ",all_reads.unibases.k" + strK )
    + String( ",all_reads.unipaths.predicted_count.k" + strK )
    + String( ",filled_reads_filt.paths.k" + strK )
    + String( ",jump_reads_ec.fastb" )
    + String( ",jump_reads_ec.distribs" )
    + String( ",}\"" );
  
  String theCommand
    = String( "RunAllPathsLG" )
    + String( " K=" + strK )
    + String( " THREADS=" + ToString( NUM_THREADS ) )
    + String( " PRE=" + PRE )
    + String( " REFERENCE_NAME=" + REFERENCE_NAME )
    + String( " DATA_SUBDIR=" + DATA_SUBDIR )
    + String( " RUN=" + RUN )
    + String( " OVERWRITE=" + String( OVERWRITE ? "True" : "False" ) )
    + String( " TARGETS= " )
    + String( " TARGETS_RUN=" + targets_run );

  RunCommand( theCommand );

}

