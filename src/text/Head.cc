/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// Head.  Read standard input until string S is encountered, copy all lines up
// to including the one with S to standard output, then exit.
//
// The following does the same thing and is faster:
// awk '{print}; /S/ {exit}'

#include "MainTools.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(S);
     EndCommandArguments;

     String line;
     while(cin)
     {    getline( cin, line );
          cout << line << "\n";
          if ( line.Contains(S) ) break;    }    }
