/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "system/System.h"

/**
 * HyperToPs
 *
 * Trivial wrapper around a system call to dot -Tps, to generate a ps
 * file from any hyper* object (either basevector or KmerPath).
 *
 * HYPER: input hyper* object (full path name)
 * PSOUT: full path name of output .ps (defaulted to HYPER.ps)
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  BeginCommandArguments;
  CommandArgument_String( HYPER );
  CommandArgument_String_OrDefault( PSOUT, "" );
  EndCommandArguments;

  // File names.
  String out_file = ( PSOUT == "" ) ? HYPER + ".ps" : PSOUT;

  // Early exit if HYPER does not exist.
  if ( ! IsRegularFile( HYPER ) ) {
    cout << Date( ) << ": input hyper not found. Exit.\n" << endl;
    return 0;
  }

  // Early exit if dot is not found.
  String comm = "dot -V > /dev/null 2>&1";
  if ( System( comm ) != 0 ) {
    cout << Date( ) << ": dot not found: run \"use graphviz\" first.\n" << endl;
    return 0;
  }
 
  // Run dot.
  cout << Date( ) << ": generating ps" << endl;
  comm =  "dot -Tps -o " + out_file + " " + HYPER;
  SystemSucceedQuiet( comm );

  // Done.
  cout << Date( ) << ": done\n" << endl;

}
