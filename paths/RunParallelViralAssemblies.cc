///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
// MakeDepend: dependency AssembleViral
// MakeDepend: dependency SelectRandomPairs

#include <omp.h>
#include "CoreTools.h"
#include "MainTools.h"
#include "TokenizeString.h"
#include "VecUtilities.h"
#include "util/RunCommand.h"

/**
 * RunParallelViralAssemblies
 *
 * Load a list file with specifics about HIV assemblies (see below),
 * and run parallel instances of AssembleViral.
 *
 * The LIST_FILE ignores empty lines, or lines starting with "#". All
 * other lines must contain exactly four comma separated Strings, as
 * in: "G15671 A07N4 1 Solexa-75632.tagged_938". In general, these
 * are: {descriptor, flowcell, lane, library). From these, an assembly
 * tree is built as follows:
 *   PRE:  <PRE>
 *   DATA: <DATA_BASE>/<flowcell>.<lane>.<library>
 *   RUN:  either <COV>x.k<K> or all.k<K>
 *
 * REF_HEAD: reference
 * GSIZE: needed to compute coverage
 * COV: downsample up to COVx coverage (or use "all")
 * FLOWCELL: if not empty, only select projects from this flowcell
 * ARGS: extra args to be passed to AssembleViral
 * NUM_THREADS: 0 means use all
 * PREPARE_ONLY: just prepare run dir
 * NOGO: just show which assemblies would be run (dry run)
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_Int( K );
  CommandArgument_String( PRE );
  CommandArgument_String( DATA_BASE );
  CommandArgument_String( REF_HEAD );
  CommandArgument_String( LIST_FILE );
  CommandArgument_Int_OrDefault( GSIZE, 10000 );
  CommandArgument_String_OrDefault( COV, "all" );
  CommandArgument_String_OrDefault( FLOWCELL, "" );
  CommandArgument_String_OrDefault( ARGS, "" );
  CommandArgument_UnsignedInt_OrDefault( NUM_THREADS, 4 );
  CommandArgument_Bool_OrDefault( PREPARE_ONLY, False );
  CommandArgument_Bool_OrDefault( NOGO, False );
  EndCommandArguments;

  // Thread control.
  NUM_THREADS = configNumThreads(NUM_THREADS);
  omp_set_num_threads( NUM_THREADS );

  // Check arguments.
  if ( COV != "all" ) {
    if ( ! COV.IsDouble( ) ) {
      cout << "Fatal error, invalid COV. Exit.\n" << endl;
      return 1;
    }
  }

  // Generate list of assemblies to run.
  vec< vec<String> > to_run;
  {
    vec<String> info;
    ifstream in( LIST_FILE.c_str( ) );
    while ( in ) {
      String line;
      getline( in, line );
      if ( ! in ) break;
      if ( line.empty( ) ) continue;
      if ( line.Contains( "#", 0 ) )continue;
      info.push_back( line );
    }

    for (int ii=0; ii<info.isize( ); ii++) {
      vec<String> tokens;
      Tokenize( info[ii], tokens );
      if ( tokens.size( ) != 4 ) continue;

      String descriptor = tokens[0];
      String flowcell = tokens[1];
      String lane = tokens[2];
      String library = tokens[3];

      if ( FLOWCELL != "" && FLOWCELL != flowcell ) continue;

      to_run.push_back( MkVec( descriptor, flowcell, lane, library ) );
    }

    sort( to_run.begin( ), to_run.end( ) );
    to_run.erase( unique( to_run.begin( ), to_run.end( ) ), to_run.end( ) );
  }
  
  // Required MB (for SelectRandomPairs).
  float mb_needed = -1.0;
  if ( COV.IsDouble( ) ) mb_needed = COV.Double( ) * (float)GSIZE / 1000000.0;

  // Dry run.
  if ( NOGO ) {
    cout << "Assemblies that would be run with NOGO=False:\n" << endl;
    for (int ii=0; ii<to_run.isize( ); ii++) {
      String descr = to_run[ii][0];
      String data = to_run[ii][1] + "." + to_run[ii][2] + "." + to_run[ii][3];
      cout << descr << "\t" <<  data << "\n";
    }
    cout << endl;
    
    cout << Date( ) << ": RunParallelViralAssemblies done" << endl;
    return 0;
  }
  
  // Loop over all listed assemblies.
  #pragma omp parallel for
  for (int ii=0; ii<to_run.isize( ); ii++) {
    String flowcell = to_run[ii][1];
    String lane = to_run[ii][2];
    String library = to_run[ii][3];
    String fll = flowcell + "." + lane + "." + library;

    String data = DATA_BASE + "/" + fll;
    String run = COV + ( COV == "all" ? "" : "x" ) + ".k" + ToString( K );
    String full_data = PRE + "/" + data;
    String full_run = full_data + "/" + run;

    String select_log = full_run + "/SelectRandomPairs.log";
    String assemble_log = full_data + "/" + run + "/AssembleViral.log";
    
    Mkpath( full_run );

    String message = Date( ) + ": starting " + fll + "/" + run + "\n";
    cout << message << flush;

    if ( COV == "all" ) {
      SymlinkForce( "../all_reads.fastb", full_run + "/reads.fastb" );
      SymlinkForce( "../all_reads.qualb", full_run + "/reads.qualb" );
      SymlinkForce( "../all_reads.pairs", full_run + "/reads.pairs" );
      SymlinkForce( "../all_reads.qltout", full_run + "/reads.qltout" );
    }
    else {
      String theCommand
	= "SelectRandomPairs READS_IN=" + full_data + "/all_reads"
	+ " READS_OUT=" + full_run + "/reads"
	+ " MB_TOTAL=" + ToString( mb_needed, 1 );
      RunCommandWithLog( theCommand, select_log );
    }

    SymlinkForce( "../params.txt", full_run + "/params.txt" );

    if ( ! PREPARE_ONLY ) {
      String theCommand
	= "AssembleViral K=" + ToString( K )
	+ " REF_HEAD=" + REF_HEAD
	+ " RUN_DIR=" + full_run
	+ ( ARGS != "" ? " " + ARGS : "" );
      RunCommandWithLog( theCommand, assemble_log );
    }

    message = Date( ) + ": done with " + fll + "/" + run + "\n";
    cout << message << flush;
  }
  
  // Done.
  cout << Date( ) << ": RunParallelViralAssemblies done" << endl;
  
}
