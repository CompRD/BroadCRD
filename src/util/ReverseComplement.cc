/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Basevector.h"

/**
 * ReverseComplement
 *
 * Reverse complement the given string.
 *
 * BASES: a string of bases
 */
int main( int argc, char *argv[] )
{
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_String( BASES );
  EndCommandArguments;

  bvec bases;
  bases.SetFromStringWithNs( BASES );
  bases.Print( cout );
  bases.ReverseComplement( );
  bases.Print( cout );
  cout << endl;

}
