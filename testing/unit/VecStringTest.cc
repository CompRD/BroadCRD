// Copyright (c) 2004 The Broad Institute

/**VecStringTest. Code to test functionality of the class vecString.
   \file VecStringTest.cc
   Some speed tests as well as verifying file read and write.

*/

#include "MainTools.h"
#include "STLExtensions.h"
#include "TaskTimer.h"
#include "VecString.h"
#include "testing/TestHelpers.h"

int main( int argc, char **argv )
{
  RunTime();
  BeginCommandArguments;
  CommandArgument_String_OrDefault( FILE, "none" );
  EndCommandArguments;

  if (FILE == "none") {
    FILE = "VecStringTestData.txt";
    vec<String> v;
    String s;
    s.resize(200);
    for (int i = 0; i != 10000; ++i) {
      StringSetRandom(s,200);
      v.push_back(s);
    }
    ofstream out(FILE.c_str());
    out << v << endl;
  }
    

  TaskTimer oldload;
  oldload.SetTrackMemory( true );
  oldload.Start();
  READ( FILE, vec<String>, oldVecString );
  oldload.Stop();
  PRINT( oldload );

  vecString newVecString;
  newVecString.reserve( oldVecString.size() );

  for ( unsigned int i = 0; i < oldVecString.size(); ++i )
    newVecString.push_back( oldVecString[i] );

  ForceAssertEq( oldVecString.size(), newVecString.size() );

  for ( unsigned int i = 0; i < oldVecString.size(); ++i )
    ForceAssertEq( oldVecString[i], newVecString[i] );

  for ( unsigned int i = 0; i < oldVecString.size(); ++i )
    oldVecString[i] = newVecString[i];

  temp_file temp( "tmp_XXXXXX" );

  newVecString.WriteAll( temp );
  
  vecString newVecStringCopy;

  TaskTimer newload;
  newload.SetTrackMemory( true );
  newload.Start();
  newVecStringCopy.ReadAll( temp );
  newload.Stop();
  PRINT( newload );

  ForceAssertEq( newVecString.size(), newVecStringCopy.size() );

  for ( size_t i = 0; i < newVecString.size(); ++i )
  {
    //ForceAssert( oldVecString[i].SelfOwned() );
    //ForceAssert( ! newVecString[i].SelfOwned() );
    //ForceAssert( ! newVecStringCopy[i].SelfOwned() );
    ForceAssertEq( newVecString[i], newVecStringCopy[i] );
  }

  //for ( int i = 0; i < newVecString.size(); ++i )
  //  ForceAssert( ! newVecString[i].SelfOwned() );
  
  TaskTimer oldsort;
  oldsort.SetTrackMemory( true );
  oldsort.Start();
  sort( oldVecString.begin(), oldVecString.end() );
  oldsort.Stop();
  PRINT( oldsort );

  TaskTimer newsort;
  newsort.SetTrackMemory( true );
  newsort.Start();
  newVecString.Sort();
  newsort.Stop();
  PRINT( newsort );

  ForceAssert( is_sorted( newVecString.begin(), newVecString.end() ) );

  for ( unsigned int i = 0; i < oldVecString.size(); ++i )
    ForceAssertEq( oldVecString[i], newVecString[i] );

  //for ( int i = 0; i < newVecString.size(); ++i )
  //  if ( newVecString[i].SelfOwned() )
  //    PRINT3( i, newVecString[i], (int)newVecString[i].SelfOwned() );
  
  //for ( int i = 0; i < newVecString.size(); ++i )
  //  ForceAssert( ! newVecString[i].SelfOwned() );
  
  cout << "All tests passed." << endl;
  
  return 0;
}
