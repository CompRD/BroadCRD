///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/**
 * Program: ChangeLibraryStats
 */
#include "MainTools.h"
#include "PairsManager.h"
#include "TokenizeString.h"
#include "VecUtilities.h"
#include <map>

int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String_Doc(HEAD_IN, "looks for HEAD_IN.pairs");
  CommandArgument_String_Doc(HEAD_OUT, "saves HEAD_OUT.pairs");
  CommandArgument_String_Doc(DATA_IN, 
    "full path name to file containing library stats (each line contains\n"
    "space separated info for one library: library_name sep sd)");
  EndCommandArguments;
  
  String pairs_file = HEAD_IN + ".pairs";
  
  cout << Date() << ": loading pairs" << endl;
  PairsManager pairs;
  if ( IsRegularFile( pairs_file ) ) pairs.Read( pairs_file );
  else FatalErr(" ERROR: did not find pairs file ");
  int nLibraries = pairs.nLibraries();
 
  Ifstream( in, DATA_IN );
  String line;
  while ( getline(in,line) ){
    vec<String> tokens;
    Tokenize( line, tokens );
    ForceAssert( tokens.size() == 3 );
    String libName = tokens[0];
    int sep        = tokens[1].Int();
    int dev        = tokens[2].Int();
    
    int libID = pairs.libraryID( libName );
    ForceAssertLt( libID, nLibraries); ForceAssertGe( libID, 0 );
    cout << "old: " << PRINT4( libName, libID, pairs.getLibrarySep(libID), pairs.getLibrarySD(libID) );
    cout << "new: " << PRINT4( libName, libID, sep, dev );
    pairs.changeLibrarySepSd( libID, sep, dev );
  }
  in.close();

  cout << Date() << " Writing pairs" << endl;
  pairs.Write( HEAD_OUT + ".pairs" );
  cout << Date() << " Done!" << endl;
}
