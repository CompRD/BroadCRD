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
#include "VecUtilities.h"
#include "feudal/BinaryStream.h"
#include "paths/TranslateAligns.h"
#include "paths/reporting/PhysicalCoverageFromAlignsCore.h"

/**
 * PhysicalCoverageFromAligns
 *
 * Compute overall physical coverage of a set of reads onto the
 * specified supers. See PhysicalCoverageFromAlignsCore.h for details.
 *
 * READS: it loads <READS>.pairs
 * ALIGNS: it loads <ALIGNS>.{qltoutlet,index}
 * SUPERS: it loads <SUPERS>.superb
 * MAX_STRETCH: used to define valid inserts
 * RADIUS: see above
 * VERBOSE: print also coverage for each super
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( READS );
  CommandArgument_String( ALIGNS );
  CommandArgument_String( SUPERS );
  CommandArgument_Double_OrDefault( MAX_STRETCH, 5.0 );
  CommandArgument_Int_OrDefault( RADIUS, 40000 );
  CommandArgument_Bool_OrDefault( VERBOSE, True );
  EndCommandArguments;

  // File names.
  String pairs_file = READS + ".pairs";
  String aligns_file = ALIGNS + ".qltoutlet";
  String index_file = ALIGNS + ".qltoutlet.index";
  String supers_file = SUPERS + ".superb";
  
  // Load.
  cout << Date( ) << ": loading pairs" << endl;
  PairsManager pairs( pairs_file );

  cout << Date( ) << ": loading aligns" << endl;
  vec<alignlet> aligns;
  BinaryReader::readFile( aligns_file, &aligns );
  
  cout << Date( ) << ": loading index" << endl;
  vec<int> index;
  BinaryReader::readFile( index_file, &index );
  
  cout << Date( ) << ": loading supers" << endl;
  vec<superb> supers;
  ReadSuperbs( supers_file, supers );

  cout << Date( ) << ": done loading\n" << endl;
  
  // Run PhysicalCoverageFromAlignsCore.
  PhysicalCoverageFromAlignsCore( cout, aligns, index, supers, pairs,
				  MAX_STRETCH, RADIUS, VERBOSE );
  
  // Done.
  cout << Date( ) << ": done" << endl;

}
