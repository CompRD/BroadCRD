/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "ParseSet.h"
#include "STLExtensions.h"

/**
 * SelectFastbSym
 *
 * Select specified entries from a fastb file, and save in ouput a
 * fastb (symmetric, ie parallel to the original) contaning exactly
 * the same number of entries as the input file, in which the only
 * non-zero entries are the ones selected.
 *
 * FASTB_IN: full path name of input fastb
 * FASTB_OUT: full path name of output fastb
 * SELECT_IDS: parsed with ParseIntSet
 */
int main( int argc, char *argv[] )
{
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_String( FASTB_IN );
  CommandArgument_String( FASTB_OUT );
  CommandArgument_String( SELECT_IDS );
  EndCommandArguments;

  vec<int> select;
  ParseIntSet( SELECT_IDS, select );
  
  cout << Date( ) << ": loading fastb" << endl;
  vecbvec bases( FASTB_IN );
  
  cout << Date( ) << ": selecting entries" << endl;
  vec<Bool> keeper( bases.size( ), False );
  for (uint ii=0; ii<select.size( ); ii++)
    keeper[ select[ii] ] = True;
  
  for (size_t ii=0; ii<bases.size( ); ii++) {
    if ( keeper[ii] ) continue;
    bases[ii].resize( 0 );
  }

  cout << Date( ) << ": saving fastb" << endl;
  bases.WriteAll( FASTB_OUT );

  cout << Date( ) << ": done" << endl;

}
