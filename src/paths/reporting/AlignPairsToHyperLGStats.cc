///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "Vec.h"
#include "feudal/BinaryStream.h"
#include "paths/Alignlet.h"
#include "paths/reporting/CLibLG.h"

/**
 * AlignPairsToHyperLGStats
 *
 * For each library, show how many pairs own reads aligning on
 * different scaffolds. It assumes AlignPairsToHyperLG has run.
 *
 * INPUT
 *   from RUN:    ALIGNS.{pairs}
 *   from SUBDIR: ALIGNS.{qltoutlet,qltoutlet.index}
 *
 * OUTPUT:
 *   sent to cout
 */
int main( int argc, char *argv[] )
{
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String_OrDefault( ALIGNS, "scaffold_reads" );
  EndCommandArguments;

  // Dir dna file names.
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;

  String pairs_file = run_dir + "/scaffold_reads.pairs";
  String aligns_file = sub_dir + "/scaffold_reads.qltoutlet";
  String index_file =  sub_dir + "/scaffold_reads.qltoutlet.index";

  // Load.
  cout << Date( ) << ": loading aligns" << endl;
  vec<alignlet> aligns;
  BinaryReader::readFile( aligns_file, &aligns );

  cout << Date( ) << ": loading index" << endl;
  vec<int> index;
  BinaryReader::readFile( index_file, &index );
  
  cout << Date( ) << ": loading pairs" << endl;
  PairsManager pairs( pairs_file );
  size_t n_pairs = pairs.nPairs( );

  // List of libs, and containers for overall stats.
  vec<String> lib_names = pairs.getLibraryNames( );
  vec<CLibLG> libs( lib_names.size( ) );
  for (int ii=0; ii<libs.isize( ); ii++) {
    libs[ii].SetLib( ii, lib_names[ii] );
    libs[ii].SetPointers( &pairs, &aligns, &index );
  }
  
  // Parse reads in pairs.
  cout << Date( ) << ": parsing " << n_pairs << " pairs" << endl;
  for (size_t pair_id=0; pair_id<n_pairs; pair_id++)
    libs[ pairs.libraryID( pair_id ) ].AddPairInfo( pair_id );

  if ( libs.size( ) < 1 ) {
    cout << "No libraries found. Exiting.\n" << endl;
    return 0;
  }
  
  // Generate and print table (start with legend).
  vec< vec<String> > table;
  table.push_back( libs[0].TableLine( True ) );
  for (size_t ii=0; ii<libs.size( ); ii++)
    table.push_back( libs[ii].TableLine( ) );

  String justif = "";
  for (int ii=0; ii<table[0].isize( ); ii++)
    justif += "r";

  cout << "\n";
  PrintTabular( cout, table, 3, justif );
  cout << "\n";
  
  // Done.
  cout << Date( ) << ": done" << endl;
  
}

