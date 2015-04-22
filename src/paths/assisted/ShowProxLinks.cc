/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "PairsManager.h"
#include "TokenizeString.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"
#include "paths/ReadLoc.h"
#include "paths/assisted/CShowLinks.h"
#include "util/RunCommand.h"

/**
 * ShowProxLinks
 *
 * Tool to show links between any two oriented contigs.
 *
 * RUN_DIR: where jump and long_jump reads are
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( CONTIGS_HEAD );
  CommandArgument_String( RUN_DIR );
  EndCommandArguments;
  
  // Dir and file names.
  String contigs_F = CONTIGS_HEAD + ".fastb";
  String locs_F = CONTIGS_HEAD + ".readlocs";
  String hits_F = CONTIGS_HEAD + ".on_reference.qlt";
  String jump_pairs_F = RUN_DIR + "/jump_reads_filt_cpd.pairs";
  String Jump_pairs_F = RUN_DIR + "/long_jump_reads_filt.pairs";

  bool use_ref = IsRegularFile( hits_F );
  bool use_jumps = IsRegularFile( jump_pairs_F );
  bool use_Jumps = IsRegularFile( Jump_pairs_F );
  if ( ! ( use_ref || use_jumps || use_Jumps ) ) {
    cout << "Fatal error: neither aligns on ref nor reads found.\n" << endl;
    return 1;
  }
  
  vec<String> needed;
  needed.push_back( contigs_F );
  needed.push_back( locs_F );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;

  // Load.
  vecbvec bases( contigs_F );

  SLibInfo ijumps;
  if ( use_jumps ) ijumps.Set( jump_pairs_F );

  SLibInfo iJumps;
  if ( use_Jumps ) iJumps.Set( Jump_pairs_F );

  vec<look_align> hits;
  if ( use_ref ) LoadLookAligns( hits_F, hits );

  String head = Basename( CONTIGS_HEAD );
  read_locs_on_disk locs_file( CONTIGS_HEAD, RUN_DIR );
  
  // Go interactive.
  CShowLinks linker( bases, ijumps, iJumps, hits, locs_file );
  linker.GoInteractive( );

  // Done.
  return 0;
  
}
