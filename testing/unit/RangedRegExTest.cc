// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#include "text/RangedRegEx.h"
#include <iostream>

int main(int argc, char** argv)
{
    using std::cout;
    using std::endl;
  const int num_tests = 13;
  const String tests[num_tests] = 
  { String( "[#1-1]" ),
    String( "paris in [#1-1] the spring" ),
    String( "paris [#1-1] in the [#1-1] spring" ),
    String( "[#1-2]" ),
    String( "[#1-9]" ),
    String( "[#1-1][#2-2]" ),
    String( "[#1-12]" ),
    String( "[#1-100]" ),
    String( "[#2345-2857]" ),
    String( "[#9-12345]" ),
    String( "[1-389]" ),
    String( "[#-10-1]" ),
    String( "[#1-2" ) } ;

  const String answers[num_tests] = 
  { String( "(1)" ),
    String( "paris in (1) the spring" ),
    String( "paris (1) in the (1) spring" ),
    String( "([12])" ),
    String( "([1-9])" ),
    String( "(1)(2)" ),
    String( "(([1-9])|(1[0-2]))" ),
    String( "(([1-9])|([1-9][0-9])|(100))" ),
    String( "((234[5-9])|(23[5-9][0-9])|(2[4-7][0-9][0-9])|(28[0-4][0-9])|(285[0-7]))" ),
    String( "((9)|([1-9][0-9])|([1-9][0-9][0-9])|([1-9][0-9][0-9][0-9])|(1[01][0-9][0-9][0-9])|(12[0-2][0-9][0-9])|(123[0-3][0-9])|(1234[0-5]))" ),
    String( "[1-389]" ),
    String( "[#-10-1]" ),
    String( "[#1-2") };

  for (int test_num = 0; test_num < num_tests; test_num++) {
    String test( tests[test_num] );
    SubstituteIntegerRanges(test);
    if ( test == answers[test_num] ) {
      cout << "Passed test " << test_num+1 << endl ;
    } else {
      cout << "Failed test " << test_num+1 << endl ;
      cout << "  '" << tests[test_num] << "' became '" 
	   << test << "'" << endl;
      cout << "  but should be '" << answers[test_num] << "'" << endl;
    }
  }

  exit(0);
}
