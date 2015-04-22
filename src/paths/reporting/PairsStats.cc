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

/**
 * PairsStats
 *
 * Minimalistic tool to print basic pairs info.
 *
 * HEAD: it loads <HEADS>.pairs
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( HEAD );
  EndCommandArguments;

  String pairs_file = HEAD + ".pairs";

  PairsManager pairs( pairs_file );

  pairs.printLibraryStats( cout );
  cout << endl;
  
}

