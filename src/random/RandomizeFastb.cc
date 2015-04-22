// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

// RandomizeFastb: take as input a given fastb file, produce as output a new
// fastb file having the same number of entries and same size entries, but having
// random bases.

#include "Basevector.h"
#include "MainTools.h"
#include "random/Random.h"

int main( int argc, char *argv[] )
{    
     RunTime( );
     
     BeginCommandArguments;
     CommandArgument_String(IN);
     CommandArgument_String(OUT);
     CommandArgument_Bool_OrDefault(CAP, False);
     EndCommandArguments;

     vecbasevector b(IN);
     for ( size_t i = 0; i < b.size( ); i++ )
     {    basevector& x = b[i];
          for ( int j = 0; j < x.isize( ); j++ )
          {    x.Set( j, randomx( ) % 4 );
               if ( CAP && j > 0 )
                    while( x[j] == x[j-1] ) x.Set( j, randomx( ) % 4 );    }    }
     b.WriteAll(OUT);    }
