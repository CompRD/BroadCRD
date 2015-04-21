/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "random/Random.h"
#include "VecUtilities.h"

#include <set>
#include <map>

int main( int argc, char** argv )
{
  srandomx(0);

  vec<long> a, b;
  set<long> a_used, b_used;
  map<long,long> truth;
  for ( int i = 0; i < 20; ++i ) {

    long newa = randomx();
    while ( a_used.count( newa ) )
      newa = randomx();
    a.push_back( newa );
    a_used.insert( newa );

    long newb = randomx();
    while ( b_used.count( newb ) )
      newb = randomx();
    b.push_back( newb );
    b_used.insert( newb );

    truth.insert( make_pair( newa, newb ) );
  }

  SortSync( a, b );

  for ( unsigned int i = 1; i < a.size(); ++i ) 
    ForceAssertLt( a[i-1], a[i] );

  for ( unsigned int i = 0; i < a.size(); ++i )
    ForceAssertEq( b[i], truth[ a[i] ] );

  SortSync( a, b, greater<long>() );

  for ( unsigned int i = 1; i < a.size(); ++i ) 
    ForceAssertGt( a[i-1], a[i] );

  for ( unsigned int i = 0; i < a.size(); ++i )
    ForceAssertEq( b[i], truth[ a[i] ] );

  cout << "Passed." << endl;
}
