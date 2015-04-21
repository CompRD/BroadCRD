// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#include "Macros.h"
#include "system/RunTime.h"

#include <iostream>
#include <getopt.h>

int main(int argc, char** argv) {
  
  RunTime();

  const int num_macros = 7;
  String macro_names[num_macros] =
  {
    String("nick"),
    String("first"),
    String("last"),
    String("middle"),
    String("middle_initial"),
    String("full"),
    String("address")
  };

  String expansions[num_macros] =
  {
    String("Jon"),
    String("${nick}athan"),
    String("Butler"),
    String("Kelly"),
    String("K."),
    String("$first $middle_initial $last"),
    String("$full, J.P. MA 02130")
  };

  const int num_tests = 15;
  String tests[num_tests] =
  {
    String("My name is $first."),
    String("My full name is $full."),
    String("$full is my full name."),
    String("My full name is ${full}."),
    String("${full} is my full name."),
    String("My address is $address."),
    String("$this is an unknown macro."),
    String("${this} is also an unknown macro."),
    String("A known macro: $full"),
    String("An unknown macro: $unknown"),
    String("A known bracketed macro: ${full}"),
    String("An unknown bracketed macro: ${unknown}"),
    String("\\$This should not be substituted."),
    String("Nor should \\$this."),
    String("My full name is $first $middle_initial $last.")
  };

  String answers[num_tests] = 
  {
    String("My name is Jonathan."),
    String("My full name is Jonathan K. Butler."),
    String("Jonathan K. Butler is my full name."),
    String("My full name is Jonathan K. Butler."),
    String("Jonathan K. Butler is my full name."),
    String("My address is Jonathan K. Butler, J.P. MA 02130."),
    String("$this is an unknown macro."),
    String("${this} is also an unknown macro."),
    String("A known macro: Jonathan K. Butler"),
    String("An unknown macro: $unknown"),
    String("A known bracketed macro: Jonathan K. Butler"),
    String("An unknown bracketed macro: ${unknown}"),
    String("$This should not be substituted."),
    String("Nor should $this."),
    String("My full name is Jonathan K. Butler.")
  };

  bool verbose, errorflag;
  verbose = errorflag = false;

  const char options[] = "v";
  char option;
  optarg = NULL; 
  while ( !errorflag && ((option = getopt(argc, argv, options)) != -1))
    switch(option) {
    case 'v':
      verbose = true;
      break;
    default:
      errorflag = true;
    }
  if ( errorflag ) {
    cerr << "Usage: " << argv[0] << " [-v]" << endl;
    exit(-1);
  }


  MacroSet macros;
  for (int macro = 0; macro < num_macros; macro++)
    macros.add(macro_names[macro], expansions[macro]);

  for (int test_num = 0; test_num < num_tests; test_num++) {
    String expanded_test;
    expanded_test = macros.expand(tests[test_num]);
    if ( expanded_test == answers[test_num] ) {
      cout << "Passed test " << test_num+1 << endl;
      if ( verbose )
	cout << "  '" << tests[test_num] << "' became '" 
	     << expanded_test << "'." << endl;
    }
    else {
      cout << "Failed test " << test_num+1 << endl;
      cout << "  '" << tests[test_num] << "' became '" << expanded_test << "'." << endl;
      cout << "  but should be '" << answers[test_num] << "'" << endl;
    }
  }

  if ( verbose ) {
    cout << "Macro listing: " << endl;
    macros.print(cout);
  }

  exit(0);
}
    
