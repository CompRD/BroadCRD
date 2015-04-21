///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "PairsManager.h"
#include "Vec.h"
#include "VecUtilities.h"
#include "feudal/BinaryStream.h"
#include "graph/Digraph.h"
#include "paths/reporting/CLinkBundle.h"
#include "paths/reporting/CSuperLinks.h"

/**
 * DumpOffsets
 *
 * Print all offsets between SUPER1 and (RC2 of) SUPER2.
 */
int main( int argc, char *argv[] )
{
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_Int( SUPER1 );
  CommandArgument_Int( SUPER2 );
  CommandArgument_Bool( RC2 );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String_OrDefault( READS, "scaffold_reads" );
  CommandArgument_String_OrDefault( ALIGNS, "scaffold_reads" );
  CommandArgument_String_OrDefault( SCAFFOLDS, "initial_scaffolds" );
  EndCommandArguments;

  // Dir and file names.
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;

  String pairs_file = run_dir + "/" + READS + ".pairs";
  String aligns_file = sub_dir + "/" + ALIGNS + ".qltoutlet";
  String index_file =  sub_dir + "/" + ALIGNS + ".qltoutlet.index";
  String supers_file = sub_dir + "/" + SCAFFOLDS + ".superb";

  // Load.
  cout << Date( ) << ": loading supers" << endl;
  vec<superb> supers;
  ReadSuperbs( supers_file, supers );
  
  cout << Date( ) << ": loading pairs" << endl;
  PairsManager pairs( pairs_file );
  size_t n_pairs = pairs.nPairs( );
  
  cout << Date( ) << ": loading aligns" << endl;
  vec<alignlet> aligns;
  BinaryReader::readFile( aligns_file, &aligns );
  
  cout << Date( ) << ": loading index" << endl;
  vec<int> index;
  BinaryReader::readFile( index_file, &index );
  
  // Find and print all offsets
  CSuperLinks linker( &pairs, &supers, &aligns, &index );
  COffset offset = linker.AllLinks( SUPER1, SUPER2, RC2 );
  offset.Print( cout, &pairs );
  
}
