///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// ColTab n1 ... nk: read from standard input, select columns n1 through nk,
// print them in order, separated by tabs.  The first column is column 1.

#include <strstream>

#include "MainTools.h"
#include "math/Functions.h"

int main( int argc, char *argv[] )
{    vec<int> cols( argc - 1 );
     for ( int i = 0; i < argc - 1; i++ )
     {    String s( argv[i+1] );
          if ( !s.IsInt( ) || s.Int( ) <= 0 )
          {    cout << "column number doesn't make sense" << endl;
               exit(1);    }
          cols[i] = s.Int( );    }
     String line;
     int n = Max(cols);
     vec<String> v(n);
     while(1)
     {    getline( cin, line );
          if ( !cin ) break;
          istrstream iline( line.c_str( ) );
          for ( int i = 0; i < n; i++ )
               v[i].resize(0);
          for ( int i = 0; i < n; i++ )
          {    if ( !iline ) break;
               iline >> v[i];
               if ( !iline ) break;    }
          for ( int j = 0; j < cols.isize( ); j++ )
               cout << ( j > 0 ? "\t" : "" ) << v[ cols[j] - 1 ];
          cout << "\n";    }    }
