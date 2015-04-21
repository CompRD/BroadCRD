/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


/// \file FeudalSize.cc: report the number of objects in a feudal file.

#include "MainTools.h"

int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String_Abbr(FILE, F);
  EndCommandArguments;

  ForceAssert( IsRegularFile(FILE) );
  cout << MastervecFileObjectCount(FILE) << " inner-vecs" << endl;
  cout << MastervecFileRawCount(FILE) << " objects in inner-vecs " << endl;

  return 0;
}
