///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Basevector.h"
#include "ParseSet.h"

/**
 * SelectSparseFastb
 *
 * Yet another module to select a subset of entries from a given
 * fastb.  Reads are loaded with SparseRead( ), and saved to a file
 * that has the exact number of entries as the input file (where only
 * the selected entries have been populated).
 *
 * IDS: parsed with ParseIntSet
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( FASTB_IN );
  CommandArgument_String( FASTB_OUT );
  CommandArgument_String( IDS );
  EndCommandArguments;

  vec<int> select;
  ParseIntSet( IDS, select );
  select.erase( unique( select.begin( ), select.end( ) ), select.end( ) );
  
  cout << Date( ) << ": loading fastb" << endl;
  vecbvec bases;
  bases.SparseRead( FASTB_IN, select, 0 );

  cout << Date( ) << ": saving fastb" << endl;
  bases.WriteAll( FASTB_OUT );

  cout << Date( ) << ": done" << endl;

}
