///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "util/RunCommand.h"
// MakeDepend: dependency BuildBackboneScaffolds
// MakeDepend: dependency BuildIdealBackbones
// MakeDepend: dependency BuildProxDigraph
// MakeDepend: dependency EvalScaffolds
// MakeDepend: dependency FillBackboneGaps
// MakeDepend: dependency InsertsCloudStats
// MakeDepend: dependency MapECJumps
// MakeDepend: dependency MapECJumpsOnUnibases
// MakeDepend: dependency MakeRcDb
// MakeDepend: dependency PrepareUnibases
// MakeDepend: dependency RunAlignlets2ReadLocs
// MakeDepend: dependency UnifiedGlobalCloud

/**
 * AssistAllPathsLG
 *
 * Assisted assembly pipeline (run SetupAssistedAssembly to generate
 * the data needed to run AssistAllPathsLG).
 */
int main( int argc, char *argv[] )
{
  RunTime( );
 
  BeginCommandArguments;

  // Standard AllPathsLG args.
  CommandArgument_String( PRE );
  CommandArgument_String( REFERENCE_NAME );
  CommandArgument_String( DATA_SUBDIR );
  CommandArgument_String( RUN );
  CommandArgument_String( SUBDIR );
  CommandArgument_Int_OrDefault( K, 96 );

  // Relative to RUN (the input reads).
  CommandArgument_String_OrDefault( UNIPATHS, "all_reads" );
  CommandArgument_String_OrDefault( FRAGS, "filled_reads_filt" );
  CommandArgument_String_OrDefault( JUMPS, "jump_reads_ec");

  // Head names of references (for assisting or evaluation), full path names.
  CommandArgument_String_OrDefault( ASSIST_REF, "" );
  CommandArgument_String_OrDefault( EVAL_REF, "" );
  
  // Reporting options.
  CommandArgument_Bool_OrDefault( INSERTS_CLOUD_STATS, True );
  
  // Create an "ideal" assembly, based on EVAL_REF (which must be given).
  CommandArgument_Bool_OrDefault( IDEAL, False );
  
  // BuildProxDigraph extra args.
  // CommandArgument_Int_OrDefault( BPD_MAX_CN, 9999999 );
  
  // Dry run, and optional start-stop args.
  CommandArgument_String_OrDefault( START, "" );
  CommandArgument_String_OrDefault( STOP, "" );
  CommandArgument_Bool_OrDefault( NOGO, False );

  EndCommandArguments;
  
  // Dir and file names.
  String strK = ToString( K );
  
  String pd = "PRE=" + PRE + " DATA=" + REFERENCE_NAME + "/" + DATA_SUBDIR;
  String pdr = pd + " RUN=" + RUN;
  String pdrs = pdr + " SUBDIR=" + SUBDIR;

  String run_dir = PRE + "/" + REFERENCE_NAME + "/" + DATA_SUBDIR + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  
  String target = IDEAL ? "ideal" : "standard";
  String full_target = sub_dir + "/" + target;
  
  // Check arguments.
  if ( IDEAL && EVAL_REF == "" ) {
    cout << "Fatal error: EVAL_REF must be provided, if IDEAL=True\n" << endl;
    return 1;
  }
  
  // List of modules to run.
  vec<String> commands;

  ///////////////////////////////////////////
  //       Starting to add commands.       //
  ///////////////////////////////////////////

  // Generate unified global cloud, and run MakeRcDb.
  String theCommand
    = "MakeRcDb"
    + String( " K=" + strK )
    + String( " " + pdr )
    + String( " READS=" + UNIPATHS )
    + String( " PATHS=unipaths" );
  commands.push_back( theCommand );

  theCommand
    = "UnifiedGlobalCloud"
    + String( " UNIPATHS=" + UNIPATHS )
    + String( " " + pdr );
  commands.push_back( theCommand );
  
  theCommand
    = "MakeRcDb"
    + String( " K=" + strK )
    + String( " " + pdr )
    + String( " READS=global_cloud" );
  commands.push_back( theCommand );
  
  String assembly_in;
  String assembly_out;
  String full_assembly_in;
  String full_assembly_out;

  assembly_out = "select";
  full_assembly_out = full_target + "/" + assembly_out;

  theCommand
    = "PrepareUnibases"
    + String( " HEAD_UNI=" + run_dir + "/" + UNIPATHS )
    + String( " OUT_DIR=" + sub_dir )
    + String( " HEAD_REF=" + ( IDEAL ? EVAL_REF : ASSIST_REF ) )
    + String( " OUT_BASE=" + target + "/" + assembly_out );
  commands.push_back( theCommand );

  theCommand
    = "MapECJumpsOnUnibases"
    + String( " RUN_DIR=" + run_dir )
    + String( " READS=" + UNIPATHS )
    + String( " JUMPS=" + JUMPS );
  commands.push_back( theCommand );
  
  theCommand
    = "BuildProxDigraph"
    + String( " HEAD_UNI=" + run_dir + "/" + UNIPATHS )
    + String( " CONTIGS_HEAD=" + full_assembly_out )
    + String( " RUN_DIR=" + run_dir )
    + String( " JUMPS=" + JUMPS );
  //+ String( " MAX_CN=" + ToString( BPD_MAX_CN ) );
  commands.push_back( theCommand );
  
  assembly_in = assembly_out;
  assembly_out = "backbones";
  full_assembly_in = full_target + "/" + assembly_in;
  full_assembly_out = full_target + "/" + assembly_out;

  theCommand
    = ( IDEAL ? "BuildIdealBackbones" : "BuildBackboneScaffolds" )
    + String( " HEAD_IN=" + full_assembly_in )
    + String( " HEAD_OUT=" + full_assembly_out );
  commands.push_back( theCommand );
  
  theCommand
    = "MapECJumps"
    + String( " CONTIGS_HEAD=" + full_assembly_out + ".contigs" )
    + String( " RUN_DIR=" + run_dir )
    + String( " JUMPS=" + JUMPS );
  commands.push_back( theCommand );
  
  theCommand
    = "RunAlignlets2ReadLocs"
    + String( " CONTIGS_HEAD=" + full_assembly_out + ".contigs" )
    + String( " RUN_DIR=" + run_dir )
    + String( " JUMPS=" + JUMPS );
  commands.push_back( theCommand );
  
  if ( INSERTS_CLOUD_STATS ) {
    theCommand
      = "InsertsCloudStats"
      + String( " " + pdrs )
      + String( " ASSEMBLY=" + target + "/" + assembly_out )
      + String( " VERBOSE=False" )
      + String( " ARCHIVE=True" );
    commands.push_back( theCommand );
  }
  
  assembly_in = assembly_out;
  assembly_out = "final";
  full_assembly_in = full_target + "/" + assembly_in;
  full_assembly_out = full_target + "/" + assembly_out;
  
  theCommand
    = "FillBackboneGaps"
    + String( " " + pdrs )
    + String( " ASSEMBLY_IN=" + target + "/" + assembly_in )
    + String( " ASSEMBLY_OUT=" + target + "/" + assembly_out );
  commands.push_back( theCommand );
  
  if ( EVAL_REF != "" ) {
    theCommand
      = "EvalScaffolds BW_ADD=3500 FORCE=True"
      + String( " LOOKUP=" + EVAL_REF + ".lookup" )
      + String( " SCAFFOLDS=" + full_assembly_out );
    commands.push_back( theCommand );
  }
  
  ///////////////////////////////////////
  //       Done adding commands.       //
  ///////////////////////////////////////
  
  // Find start and stop.
  int start = StartCommand( commands, START );
  int stop = StopCommand( commands, STOP );
  
  // Run modules.
  vec<String> short_commands;
  ShortCommands( commands, short_commands );

  if ( NOGO ) {
    vec< vec<String> > table;
    vec<String> line;
    for (int ii=start; ii<=stop; ii++) {
      line.clear( );
      line.push_back( short_commands[ii] );
      line.push_back( "run as:" );
      line.push_back( commands[ii] );
      table.push_back( line );
    }
    PrintTabular( cout, table, 2 );
    cout << endl;
    return 0;
  }
  
  for (int ii=start; ii<=stop; ii++) {
    if ( ! RunCommandBool( commands[ii] ) ) {
      String err = "Fatal error: " + short_commands[ii] + " failed"; 
      cout << "\n " << err << "\n" << endl;
      return 0;
    }
  }
  cout << endl;
  
  // Done.  
  cout << Date( ) << ": AssistAllPathsLG done" << endl;
  
}

