// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

// Comma: read an integer from standard input and put commas in it every three 
// digits.  This is stupid -- it only handles integers which fit in 8 bytes.

#include "MainTools.h"

int main( )
{    longlong x;
     cin >> x;
     if ( x < 0 ) 
     {    cout << "-";
          x = -x;    }
     vec<int> digits;
     while(1)
     {    digits.push_back( x % 10 );
          if ( x < 10 ) break;
          x /= 10;    }
     for ( int i = digits.isize( ) - 1; i >= 0; i-- )
     {    cout << digits[i];
          if ( i > 0 && (i % 3) == 0 ) cout << ",";    }
     cout << "\n";    }
