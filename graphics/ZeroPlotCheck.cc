/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// This just checks if the input to PlotPoints has all zeros in the Y axis.
//       *** Returns a 1 if the plot has all zeros. ***
// This is meant to be used in, for example, a Perl script as a special case
// check.  Use with NH=True.  Works on both double and integer input files.

#include "MainTools.h"
#include <fstream>

int main( int argc, char *argv[] ) 
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(HEAD);
     EndCommandArguments;

     std::ifstream is(HEAD.c_str());
     int all0 = 1;
     double x, y;
     while ( (is >> x >> y) )
     {
         if ( y != 0 )
         {
             all0 = 0;
             break;
         }
     }
     ForceAssert(is.eof());
     std::cout << all0 << std::endl;
     return all0;
}
