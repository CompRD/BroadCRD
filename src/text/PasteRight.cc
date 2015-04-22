// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

// PasteRight: for each line of FILE1, pad with blanks on right to size N.
// Then paste the lines of FILE2 on the right.

#include "MainTools.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(FILE1);
     CommandArgument_String(FILE2);
     CommandArgument_UnsignedInt(N);
     EndCommandArguments;

     Ifstream( in1, FILE1 );
     Ifstream( in2, FILE2 );

     String line1, line2;
     while(1)
     {    getline( in1, line1 );
          getline( in2, line2 );
          if ( !in1 || !in2 ) break;
          cout << line1;
          for ( int i = 0; i < (int) N - (int) line1.size( ); i++ )
               cout << " ";
          cout << line2 << "\n";    }    }
