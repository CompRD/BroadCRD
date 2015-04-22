/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// SortParagraphs: read paragraphs from input file, sort lexicographically, print.
//
// If UNIQUE=True, remove duplicated paragraphs and print paragraph multiplities.

#include "FastIfstream.h"
#include "MainTools.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(INPUT);
     CommandArgument_String(OUTPUT);
     CommandArgument_Bool_OrDefault(UNIQUE, False);
     EndCommandArguments;

     fast_ifstream in(INPUT);
     Ofstream( out, OUTPUT );

     vec<String> paragraphs;
     ReadParagraphs( in, paragraphs );

     Sort(paragraphs);
     for ( int i = 0; i < paragraphs.isize( ); i++ )
     {    if ( !UNIQUE ) out << paragraphs[i] << "\n";
          else
          {    int j;
               for ( j = i + 1; j < paragraphs.isize( ); j++ )
                    if ( paragraphs[j] != paragraphs[i] ) break;
               out << "COUNT = " << j - i << "\n";
               out << paragraphs[i] << "\n";
               i = j - 1;    }
          if ( i < paragraphs.isize( ) - 1 ) out << "\n";    }    }
