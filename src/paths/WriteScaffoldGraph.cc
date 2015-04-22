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
#include "PairsManager.h"
#include "Superb.h"
#include "Vec.h"
#include "paths/Alignlet.h"
#include "paths/Sepdev.h"
#include "paths/BuildScaffoldGraph.h"
#include "paths/reporting/CLinkBundle.h"
#include "util/RunCommand.h"
#include "feudal/BinaryStream.h"
// MakeDepend: dependency ScaffoldGraphToGnuplot

/**
 * WriteScaffoldGraph
 *
 * Run BuildScaffoldGraph and save the generated scaffold graph to
 * file.
 *
 * GNUPLOT_READY: if true, generate also output for links-plot.pl
 */
int main( int argc, char *argv[] )
{
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String_OrDefault( READS, "scaffold_reads" );
  CommandArgument_String_OrDefault( ALIGNS, "scaffold_reads" );
  CommandArgument_String_OrDefault( SCAFFOLDS, "initial_scaffolds.superb" );
  CommandArgument_String_OrDefault( OUT_GRAPH, "initial_scaffolds.graph" );
  CommandArgument_Bool_OrDefault( GNUPLOT_READY, True );
  EndCommandArguments;

  // Dir and file names.
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;

  String pairs_file = run_dir + "/" + READS + ".pairs";
  String aligns_file = sub_dir + "/" + ALIGNS + ".qltoutlet";
  String index_file =  sub_dir + "/" + ALIGNS + ".qltoutlet.index";
  String supers_file = sub_dir + "/" + SCAFFOLDS;
  String graph_file = sub_dir + "/" + OUT_GRAPH;
  String head_file = graph_file.SafeBefore( ".graph" );

  // Load.
  cout << Date( ) << ": loading pairs" << endl;
  PairsManager pairs( pairs_file );
  size_t n_pairs = pairs.nPairs( );

  cout << Date( ) << ": loading supers" << endl;
  vec<superb> supers;
  ReadSuperbs( supers_file, supers );
  
  cout << Date( ) << ": loading aligns" << endl;
  vec<alignlet> aligns;
  BinaryReader::readFile( aligns_file, &aligns );

  cout << Date( ) << ": loading index" << endl;
  vec<int> index;
  BinaryReader::readFile( index_file, &index );
  
  // Generate initial graph.
  cout << Date( ) << ": generating scaffolds graph" << endl;
  digraphE<sepdev> unused;
  digraphE<CLinkBundle> graph;
  BuildScaffoldGraph( pairs, supers, aligns, index, unused, &graph );
  
  // Save.
  cout << Date( ) << ": saving graph" << endl;
  BinaryWriter::writeFile( graph_file, graph );
  
  if ( GNUPLOT_READY ) {
    cout << Date( ) << ": saving graph.txt" << endl;
    String theCommand = "ScaffoldGraphToGnuplot HEAD=" + head_file;
    RunCommandWithLog( theCommand, "/dev/null" );
  }

  // Done.
  cout << Date( ) << ": done" << endl;

}
