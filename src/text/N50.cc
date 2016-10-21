///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Read from standard input and return N50.

#include <strstream>

#include "MainTools.h"
#include "math/Functions.h"

int main( int argc, char *argv[] )
{    vec<int> all;
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
     {    cout << "No entries, can't compute N50" << endl;
          exit(1);    }
     cout << N50(all) << endl;    }

