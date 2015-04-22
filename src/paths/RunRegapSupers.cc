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
#include "VecUtilities.h"
#include "feudal/BinaryStream.h"
#include "math/Functions.h"
#include "paths/Alignlet.h"
#include "paths/FixScaffoldsCore.h"
#include "paths/OffsetDistribution.h"
#include "paths/ReadLoc.h"
#include "paths/RegapSupers.h"
#include "paths/SaveScaffoldGraph.h"
#include "paths/ScaffoldsUtils.h"

/**
 * RunRegapSupers
 *
 * Run RegapSupers on a specified assembly (it may require the library
 * distributions).  No output is generated - this is testing code.
 *
 * INPUT:
 *   <READS>.{distribs,pairs}
 *   <ALIGNS>.qltoutlet.{,index}
 *   <SUPERS>.superb
 *
 * OUTPUT:
 *   log (sent to cout)
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( READS );
  CommandArgument_String( ALIGNS );
  CommandArgument_String( SUPERS );
  CommandArgument_Bool_OrDefault( VERBOSE, False );
  CommandArgument_Bool_OrDefault( LIB_DIST, True );
  CommandArgument_Int_OrDefault( MAX_OVERLAP, 15000 );
  EndCommandArguments;

  String distribs_file = READS + ".distribs";
  String pairs_file = READS + ".pairs";
  String aligns_file = ALIGNS + ".qltoutlet";
  String index_file =  ALIGNS + ".qltoutlet.index";
  String supers_file = SUPERS + ".superb";
  
  vec<IntDistribution> distribs;
  if ( LIB_DIST ) {
    cout << Date( ) << ": loading library distributions" << endl;
    size_t nlibs = 0;
    BinaryReader reader( distribs_file.c_str( ) );
    reader.read( &nlibs );
    distribs.resize( nlibs );
    reader.readItr( distribs.begin( ), distribs.end( ) );
  }

  cout << Date( ) << ": loading supers" << endl;
  vec<superb> supers;
  ReadSuperbs( supers_file, supers );
  
  cout << Date( ) << ": loading aligns" << endl;
  vec<alignlet> aligns;
  BinaryReader::readFile( aligns_file, &aligns );
  
  cout << Date( ) << ": loading index" << endl;
  vec<int> index;
  BinaryReader::readFile( index_file, &index );
  
  cout << Date( ) << ": loading pairs" << endl;
  PairsManager pairs( pairs_file );
  size_t n_pairs = pairs.nPairs( );

  cout << Date( ) << ": done loading\n" << endl;
  
  RegapSupers( cout, supers, pairs, aligns, index, VERBOSE, MAX_OVERLAP,
	       1, 12000, True, ( LIB_DIST ? &distribs : 0 ) );
  
}

