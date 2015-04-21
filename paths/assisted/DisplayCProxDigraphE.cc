/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "graph/Digraph.h"
#include "paths/assisted/CProx.h"
#include "paths/assisted/CDisplayProxDigraph.h"
#include "paths/assisted/Proxgraph.h"
#include "util/RunCommand.h"

/**
 * DisplayCProxDigraphE
 *
 * Interactive tool to display portions of a given digraphE<Cprox>.
 *
 * Remark: dot must be in your path.
 *
 * INPUT
 *   <HEAD>{<CONTIGS>.fastb,.CProx.DigraphE}
 */
int main( int argc, char *argv[] )
{
  
  BeginCommandArguments;
  CommandArgument_String( HEAD );
  CommandArgument_String_OrDefault( CONTIGS, "" );
  EndCommandArguments;

  String bases_file = HEAD + CONTIGS + ".fastb";
  String graph_file = HEAD + ".CProx.digraphE";
    
  vec<String> needed;
  needed.push_back( bases_file );
  needed.push_back( graph_file );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;
  
  if ( StringOfOutput( "which dot" ) == "no" ) {
    cout << "Fatal error: dot is not installed on this system." << endl;
    return 0;
  }
  
  cout << Date( ) << ": loading bases" << endl;
  vecbvec bases( bases_file );

  cout << Date( ) << ": loading graph" << endl;
  proxgraph graph;
  BinaryReader::readFile( graph_file, &graph );

  cout << Date( ) << ": generating graph maps" << endl;
  CDisplayProxDigraph displayer( graph, bases );

  cout << Date( ) << ": ready - type h for help\n" << endl;
  displayer.GoDisplay( );

  cout << Date( ) << ": done" << endl;

}
