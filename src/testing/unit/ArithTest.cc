// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

// ArithTest.  If "Float" as defined in Arith.h works correctly, this code should
// yield the same answer on alpha, ia64, and i686 platforms.  The converse is NOT
// the case: this is in no sense a complete test of floating point arithmetic
// behavior.
//
// You could change Float to Double, and try to define a Double class that also
// satisfies this test, but I have not been able to do so.

#include "math/Arith.h"
#include "MainTools.h"
#include "random/Random.h"

int main( ) 
{   int answer = 0;
    for ( int j = 0; j < 1000000; j++ )
    {    Float x1 = Float(1) - Float(1)/Float( randomx( ) % 100000 + 1 );
         Float x2 = Float(1) / x1;
         for (int i = 0; i < 10; i++)
         {    if ( x1*x2 >= Float(1) ) answer += randomx( ) % 100;
              if ( answer > 1000000 ) answer /= 10;
              x1 *= x1;
              x2 *= x2;    }    }
    PRINT(answer);    } 
