///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "String.h"

/**
 * ParseEvaluateViralAssemblyLogs
 *
 * Gather stats from multiple runs of AssembleViral. See
 * RunParallelViralAssemblies for a description of LIST_FILE and
 * assumptions on the directory tree.
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( LIST_FILE );
  CommandArgument_String( ASSEMBLIES_DIR );
  CommandArgument_String( COMMON_RUN );
  CommandArgument_String_OrDefault( COMMON_LOG, "AssembleViral.log" );
  EndCommandArguments;

  // Parse LIST_FILE, generate a list of (existing) log files.
  vec<String> log_files;
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
      String fll = flowcell + "." + lane + "." + library;
      String run_dir = ASSEMBLIES_DIR + "/" + fll + "/" + COMMON_RUN;
      String log_file = run_dir + "/" + COMMON_LOG;
      if ( ! IsRegularFile ( log_file ) ) continue;

      log_files.push_back( log_file );
    }
  }


  



}
