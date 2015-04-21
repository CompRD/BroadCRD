///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Read from standard input and return median of column 1.

#include <strstream>

#include "MainTools.h"

int main( int argc, char *argv[] )
{    vec<double> all;
     String line;
     double x;
     while(1)
     {    getline( cin, line );
          if ( !cin) break;
          istrstream iline( line.c_str( ) );
          iline >> x;
          all.push_back(x);    }
     Sort(all);
     if ( all.size( ) == 0 ) 
     {    cout << "No entries, can't compute median" << endl;
          exit(1);    }
     cout << all[ all.size( ) / 2 ] << endl;    }

