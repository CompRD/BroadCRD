// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

// TracebackTest: determine if under various scenarios where
// TracebackThisProcess is called, a usable traceback results.

// CASE: 1 or 2, at present.

// Results for gcc 2.95.2 were obtained by copying
// testing/{TraceTest.cc, TraceTest2.cc, TraceTest2.h}
// to a "cvs -r pre_333 checkout" of Arachne.

#include "MainTools.h"
#include "testing/TraceTest2.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_UnsignedInt(CASE);
     EndCommandArguments;

     ForceAssert( CASE == 1 || CASE == 2 );

     if ( CASE == 1 )
     {
          vec<int> v(1);
          v[0] = 12;
          cout << v[5] << endl;
          cout << EvalVec(v, -1) << endl;

          // gcc 3.3.3 alpha results: 
          // 0. vec<int>::operator[](unsigned long) const, in stl_vector.h:415
          // 1. EvalVec(vec<int> const&, int), in TraceTest2.cc:6
          // 2. main, in ostream:193

          // gcc 3.3.3 ia64 results (O3 same as O2):
          // 0. vec<int>::operator[](unsigned long) const, in Vec.h:50
          // 1. EvalVec(vec<int> const&, int), in TraceTest2.cc:12
          // 2. main, in ostream:193

          // gcc 2.95.2 alpha results:
          // 0. vec<int>::operator[](unsigned long) const, in Vec.h:52
          // 1. EvalVec(vec<int> const &, int), in TraceTest2.cc:6
          // 2. main, in TraceTest.cc:21

               }

     if ( CASE == 2 )
     {
          vec<int> v;
          cout << EvalVec2( v, 5 ) << endl;

          // gcc 3.3.3 alpha results: 
          // 0. vec<int>::operator[](unsigned long) const, in stl_vector.h:415
          // 1. EvalVec2(vec<int> const&, int), in ReadLocation.h:103
          // 2. main, in ostream:193

          // gcc 3.3.3 ia64 results (O3 same as O2):
          // 0. vec<int>::operator[](unsigned long) const, in Vec.h:50
          // 1. EvalVec2(vec<int> const&, int), in TraceTest2.cc:12
          // 2. main, in ostream:193

          // gcc 2.95.2 alpha results:
          // 0. vec<int>::operator[](unsigned long) const, in Vec.h:52
          // 1. EvalVec2(vec<int> const &, int), in TraceTest2.cc:11
          // 2. main, in TraceTest.cc:45
     
               }    }

     
