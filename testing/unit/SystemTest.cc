// Copyright (c) 2005 Broad Institute of MIT and Harvard
//

// Test functions in System.h.

#include "MainTools.h"

bool TestDirname( const String& path, const String& expectedDirname )
{
  static int testNum = 0;
  ++testNum;

  String observedDirname = Dirname( path );

  if ( observedDirname == expectedDirname )
    return true;
  else
  {
    cout << "Test #" << testNum << " FAILED:"
         << " Dirname(" << path << ") returned '" << observedDirname << "',"
         << " but should be '" << expectedDirname << "'." << endl;
    return false;
  }
}

int main( int argc, char** argv )
{
  RunTime();

  int failures = 0;
  
  // Test Dirname().
  {
    if ( ! TestDirname( "/", "/" ) ) ++failures;
    if ( ! TestDirname( "asdf", "." ) ) ++failures;
    if ( ! TestDirname( "/asdf", "/" ) ) ++failures;
    if ( ! TestDirname( "//asdf", "/" ) ) ++failures;
    if ( ! TestDirname( "///asdf", "//" ) ) ++failures;
    if ( ! TestDirname( "/asdf/", "/" ) ) ++failures;
    if ( ! TestDirname( "/asdf/qwerty", "/asdf" ) ) ++failures;
    if ( ! TestDirname( "/asdf//qwerty", "/asdf/" ) ) ++failures;
    if ( ! TestDirname( "/asdf/qwerty/", "/asdf" ) ) ++failures;
    if ( ! TestDirname( "/asdf/qwerty//", "/asdf" ) ) ++failures;
  }

  if ( failures )
    cout << failures << " tests FAILED." << endl;
  else
    cout << "All tests PASSED." << endl;
  
  return ( failures > 0 );
}
