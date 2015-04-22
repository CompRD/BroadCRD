///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"

/**
 * FindIndex
 *
 * Find index of given int in a vec<int> saved with the macro WRITE.
 * The input file FILE does not need to be sorted.
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( FILE );
  CommandArgument_Int( OBJECT);
  EndCommandArguments;
  
  READ( FILE, vec<int>, data );
  for (size_t ii=0; ii<data.size( ); ii++) {
    if ( data[ii] == OBJECT ) {
      cout << "Index of " << OBJECT << ": " << ii << "\n" << endl;
      return 0;
    }
  }

  cout << "Index of " << OBJECT << ": not found\n" << endl;
  return 0;

}
