 // Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
 //
 //
 //
 // Log10Problem: log10 produces different results on
 // alpha/g++(2.95.2|3.3.3) and ia64/gcc3.3.3.
 //
 // On alpha/gcc2.95.2:
 // test_double = 0x465f8def8808b000
 // log_test =    0x403f000000000000
 // NO
 //
 // On alpha/gcc3.3.3:
 // test_double = 0x465f8def8808b000
 // log_test =    0x403f000000000000
 // NO
 //
 // On ia64:
 // test_double = 0x465f8def8808b000
 // log_test =    0x403effffffffffff
 // YES
 //
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;

int main( )
{
  double test_double = 9.99999999999995910349965e+30;
  double log_test = log10( test_double );

  // We print out the byte-representations of the doubles.
  printf( "test_double = %#lx\n", *((unsigned long*) &test_double) );
  printf( "log_test =    %#lx\n", *((unsigned long*) &log_test) );

  printf( "%s\n", ( log_test < 31.0 ? "YES" : "NO" ) );

  return 0;
}
