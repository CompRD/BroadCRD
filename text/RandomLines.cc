///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Extract random lines from standard input.

#include "MainTools.h"
#include "random/Random.h"

int main( int argc, char *argv[] )
{    
     RunTime( );

     BeginCommandArguments; 
     CommandArgument_Double(FRAC);
     EndCommandArguments; 

     int top = int( round( FRAC * 1000000.0 ) );
     String line;
     while(1)
     {    getline( cin, line );
          if ( !cin) break;
          if ( randomx( ) % 1000000 < top ) cout << line << "\n";    }    }
