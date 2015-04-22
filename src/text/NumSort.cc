// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

// NumSort: read from IN and sort lines using cmp_numeric, write to OUT.

#include "MainTools.h"
#include "FastIfstream.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(IN);
     CommandArgument_String(OUT);
     EndCommandArguments;

     int nlines = LineCount(IN);
     vec<String> lines;
     lines.reserve(nlines);
     fast_ifstream in(IN);
     String line;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          lines.push_back(line);    }
     sort( lines.begin( ), lines.end( ), cmp_numeric );
     Ofstream( out, OUT );
     for ( int i = 0; i < lines.isize( ); i++ )
          out << lines[i] << "\n";    }
