// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

// ToMatrix: read standard input, put in columns with SEP spaces between them and
// JUST justification (default: left), and write to standard output.

#include <strstream>

#include "MainTools.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_UnsignedInt_OrDefault(SEP, 2);
     CommandArgument_String_OrDefault(JUST, "");
     EndCommandArguments;

     vec< vec<String> > rows;
     String line;
     while(1)
     {    getline( cin, line );
          if ( !cin ) break;
          istrstream iline( line.c_str( ) );
          static vec<String> row;
          row.clear( );
          while(iline)
          {    static String s;
               iline >> s;
               row.push_back(s);    }
          rows.push_back(row);    }
     PrintTabular( cout, rows, SEP, JUST );    }
