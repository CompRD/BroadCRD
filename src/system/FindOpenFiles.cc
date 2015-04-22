/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// FileOpenFiles.  Parse output of "strace -e trace=open,close <command>"
// to find the files that were left open.  Barely works.

#include <map>

#include "FastIfstream.h"
#include "MainTools.h"
#include "Set.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(FILE);
     EndCommandArguments;

     fast_ifstream in(FILE);
     set<int> opens;
     map<int, String> openline;
     String line;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          if ( line.Contains( "open" ) && line.Contains( "=" ) )
          {    int open = line.RevAfter( " = " ).Int( );
               opens.insert(open);
               openline[open] = line;    }
          else 
          {    if ( line.Contains( "close(" ) 
                    && line.After( "close(" ).Contains( ")" ) )
               {    int close = line.Between( "close(", ")" ).Int( );
                    opens.erase(close);    }    }    }
     for ( set<int>::iterator i = opens.begin( ); i != opens.end( ); i++ )
          cout << openline[*i] << "\n";    }
