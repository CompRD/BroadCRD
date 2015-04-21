// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

// PullInterval: extract all the lines in a file between the first line containing
// string S1 and the first line after it containing string S2.

#include "MainTools.h"
#include "FastIfstream.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(FILE);
     CommandArgument_String(S1);
     CommandArgument_String(S2);
     EndCommandArguments;

     fast_ifstream in(FILE);
     String line;
     Bool found1 = False;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) )
          {    cout << "MATCH FAILED\n";
               exit(1);    }
          if ( !found1 && line.Contains(S1) )
          {    cout << line << "\n";
               found1 = True;    }
          else if (found1) 
          {    cout << line << "\n";
               if ( line.Contains(S2) ) break;    }    }    }
