// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

// DiffFiles: find all lines in the first file which are not in the
// second file.
//
// I assume that this can be done with a shorter shell script.  If you know
// how (and it's not slower), please feel free to replace this by the script.

#include <algorithm>

#include "MainTools.h"
#include "FastIfstream.h"

int main( int argc, char *argv[] )
{
     BeginCommandArguments;
     CommandArgument_String(FILE1);
     CommandArgument_String(FILE2);
     EndCommandArguments;

     vec<String> v1( LineCount(FILE1) ), v2( LineCount(FILE2) );
     fast_ifstream in1(FILE1), in2(FILE2);
     String s;
     int p1 = 0, p2 = 0;
     while(1)
     {    getline( in1, s );
          if ( in1.fail( ) ) break;
          v1[ p1++ ] = s;    }
     while(1)
     {    getline( in2, s );
          if ( in2.fail( ) ) break;
          v2[ p2++ ] = s;    }
     sort( v2.begin( ), v2.end( ) );
     for ( int i = 0; i < (int) v1.size( ); i++ )
          if ( !BinMember( v2, v1[i] ) ) cout << v1[i] << "\n";    }
