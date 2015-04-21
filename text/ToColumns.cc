// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

// ToColumns: read standard input, treating lines as "words".  Reformat the words
// so that they are in columns, with 2 spaces between columns, so that the page is
// at most 80 spaces wide.

#include "MainTools.h"
#include "math/Functions.h"

const int page_width = 80, col_sep = 2;

int main( )
{    vec<String> lines;
     String line;
     unsigned int max_width = 0;
     while(1)
     {    getline( cin, line );
          if ( !cin ) break;
          max_width = Max( max_width, line.size( ) );
          lines.push_back(line);    }
     int n = lines.size( );
     int ncols = ( page_width + col_sep ) / ( max_width + col_sep );
     if ( ncols <= 1 )
     {    for ( int i = 0; i < n; i++ )
               cout << lines[i] << "\n";
          return 0;    }
     int count = 0, nrows = (n + ncols - 1) / ncols;
     vec< vec<String> > rows( nrows, vec<String>(ncols) );
     for ( int i = 0; i < nrows; i++ )
     {    for ( int j = 0; j < ncols; j++ )
               if ( count < lines.isize( ) ) rows[i][j] = lines[count++];    }
     PrintTabular( cout, rows, 2 );    }
