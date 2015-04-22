/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"

int main( int argc, char** argv ) {
  RunTime();

  BeginCommandArguments;
  CommandArgument_String_OrDefault(TEST,"");
  CommandArgument_Bool_OrDefault(USE_STRING,False);
  EndCommandArguments;

  PRINT2( TEST, (bool)USE_STRING );

  if ( TEST == "" ) {
    cout << "Shouldn't assert at line " << __LINE__ + 1 << "." << endl;
    ForceAssert( true );
    cout << "Should assert at line " << __LINE__ + 1 << "." << endl;
    ForceAssert( false );
  } 

  else if ( TEST == "eq" ) {
    if ( USE_STRING ) {
      String A = "this";
      String B = "this";
      String C = "that";
      cout << "Shouldn't assert at line " << __LINE__ + 1 << "." << endl;
      ForceAssertEq( A, B );
      cout << "Should assert at line " << __LINE__ + 1 << "." << endl;
      ForceAssertEq( A, C );
    } else {
      longlong a = 0, b = 0, c = 1;
      cout << "Shouldn't assert at line " << __LINE__ + 1 << "." << endl;
      ForceAssertEq( a, b );
      cout << "Should assert at line " << __LINE__ + 1 << "." << endl;
      ForceAssertEq( a, c );
    }
  }

  else if ( TEST == "ne" ) {
    if ( USE_STRING ) {
      String A = "this";
      String B = "that";
      String C = "this";
      cout << "Shouldn't assert at line " << __LINE__ + 1 << "." << endl;
      ForceAssertNe( A, B );
      cout << "Should assert at line " << __LINE__ + 1 << "." << endl;
      ForceAssertNe( A, C );
    } else {
      longlong a = 0, b = 1, c = 0;
      cout << "Shouldn't assert at line " << __LINE__ + 1 << "." << endl;
      ForceAssertNe( a, b );
      cout << "Should assert at line " << __LINE__ + 1 << "." << endl;
      ForceAssertNe( a, c );
    }
  }

  else if ( TEST == "lt" ) {
    if ( USE_STRING ) {
      cout << "No string test for " << TEST << "." << endl;
      TracebackThisProcess();
    }
    longlong a = 0, b = 1, c = 0;
    cout << "Shouldn't assert at line " << __LINE__ + 1 << "." << endl;
    ForceAssertLt( a, b );
    cout << "Should assert at line " << __LINE__ + 1 << "." << endl;
    ForceAssertLt( a, c );
  }

  else if ( TEST == "gt" ) {
    if ( USE_STRING ) {
      cout << "No string test for " << TEST << "." << endl;
      TracebackThisProcess();
    }
    longlong a = 1, b = 0, c = 1;
    cout << "Shouldn't assert at line " << __LINE__ + 1 << "." << endl;
    ForceAssertGt( a, b );
    cout << "Should assert at line " << __LINE__ + 1 << "." << endl;
    ForceAssertGt( a, c );
  }

  else if ( TEST == "le" ) {
    if ( USE_STRING ) {
      cout << "No string test for " << TEST << "." << endl;
      TracebackThisProcess();
    }
    longlong a = 1, b = 2, c = 0;
    cout << "Shouldn't assert at line " << __LINE__ + 1 << "." << endl;
    ForceAssertLe( a, b );
    cout << "Should assert at line " << __LINE__ + 1 << "." << endl;
    ForceAssertLe( a, c );
  }

  else if ( TEST == "ge" ) {
    if ( USE_STRING ) {
      cout << "No string test for " << TEST << "." << endl;
      TracebackThisProcess();
    }
    longlong a = 1, b = 0, c = 2;
    cout << "Shouldn't assert at line " << __LINE__ + 1 << "." << endl;
    ForceAssertGe( a, b );
    cout << "Should assert at line " << __LINE__ + 1 << "." << endl;
    ForceAssertGe( a, c );
  }

  else {
    cout << "Unknown test " << TEST << "." << endl;
    TracebackThisProcess();
  }
}

