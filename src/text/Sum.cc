///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Read from standard input and add up column 1.

#include <strstream>

#include "MainTools.h"

int main( int argc, char *argv[] )
{    longlong sum = 0;
     int shift = 0;
     String line, x;
     while(1)
     {    getline( cin, line );
          if ( !cin) break;
          istrstream iline( line.c_str( ) );
          iline >> x;
          longlong mult = 1;
          if ( x.Contains( "k", -1 ) ) 
          {    mult = 1000;
               x = x.RevBefore( "k" );    }
          else if ( x.Contains( "m", -1 ) ) 
          {    mult = 1000000;
               x = x.RevBefore( "m" );    }
          else if ( x.Contains( "g", -1 ) ) 
          {    mult = 1000000000;
               x = x.RevBefore( "g" );    }
          int s = 0;
          if ( x.Contains( "." ) ) s = x.After( "." ).size( );
          longlong y;
          if ( x.Contains( "." ) ) y = ( x.Before( "." ) + x.After( "." ) ).Int( );
          else y = x.Int( );
          while ( shift < s )
          {    sum *= 10;
               shift++;    }
          while ( s < shift )
          {    y *= 10;
               s++;    }
          sum += mult * y;    }
     if ( shift == 0 ) cout << ToStringAddCommas(sum) << "\n";
     else
     {    String d = ToString(sum);
          cout << ToStringAddCommas( d.substr( 0, d.isize( ) - shift ).Int( ) );
          String frac = d.substr( d.isize( ) - shift, shift );
          if ( frac != "0" ) cout << "." << frac;    
          cout << "\n";    }    }
