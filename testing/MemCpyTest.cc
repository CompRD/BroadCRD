// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
//


// MemCpyTest: This code establishes that memcpy'ing 2^31 or more bytes
// is not safe under g++ (or at least under g++ 2.95.2, with whatever
// libraries we have, as of 11/6/01).  It is an open question whether
// such copying is safe under cxx.

// When compiled with g++ -O3 on an Alpha, this code fails.  However, when 
// compiled with cxx -O5, it runs fine.

// output: memcpy failed, i = 645813112, after about 30 seconds

#include <vector>
#include <iostream>
#include <string.h>
#include "system/RunTime.h"

using namespace std;

class nine_int {
     int x[9];
};

template<class T> inline void Destroy( vector<T> & v )
{    v.clear( );
     vector<T> tmp = v;
     v.swap(tmp);    }

int main( )
{
     RunTime();

     // n1 = 2^31
     long n1 = (long) 2147483 * (long) 1000 + (long) 648;

     long n2 = 775557972;

     char* x = new char[n1];
     for ( long i = 0; i < n1; i++ )
          x[i] = 0;

     vector<nine_int> locs, locs_orig;
     int unpp_size = 1512954;
     int locs_size = 16426299, locs_orig_size = 16164684;
     locs.reserve( locs_size + unpp_size );
     locs_orig.reserve( locs_orig_size + unpp_size );

     Destroy(locs);
     Destroy(locs_orig);

     char* y = new char[n1 + n2];

     memcpy( y, x, n1 );

     for ( long i = 0; i < n1; i++ )
          if ( y[i] != 0 )
          {    cout << "memcpy failed, i = " << i << endl;
               break;    }    }
