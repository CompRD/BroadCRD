/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// MakeFourControls.  Generate four N-base sequences that are random, subject to
// the constraints that:
// (a) at each position, every base occurs exactly once;
// (b) no base occurs twice in a row on a given sequence.
// Theoretically, this code could get stuck in an infinite loop.

#include "Basevector.h"
#include "MainTools.h"
#include "random/Random.h"
#include "random/Shuffle.h"

int main( int argc, char *argv[] ) 
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_Int(N);
     EndCommandArguments;

     vec<String> x(4);
     for ( int i = 0; i < N; i++ )
     {    
          // Shuffle the elements of x.

          vec<int> sh(4);
          int seed = randomx( );
          Shuffle( 4, sh, seed );
          vec<String> y(4);
          PRINT4( sh[0], sh[1], sh[2], sh[3] );
          for ( int j = 0; j < 4; j++ )
               y[j] = x[ sh[j] ];
          x = y;

          // Append to each x.

          int base1, base2, base3, base4;
          int count = 0;
          stuck:
          if ( ++count == 100 )
          {    i -= 2;
               for ( int j = 0; j < 4; j++ )
                    x[j].resize( x[j].size( ) - 1 );
               continue;    }
          int count1 = 0, count2 = 0, count3 = 0, count4 = 0;
          while(1)
          {    base1 = randomx( ) % 4;
               if ( ++count1 == 100 ) goto stuck;
               if ( i == 0 || as_base(base1) != x[0][i-1] ) break;    }
          while(1)
          {    base2 = randomx( ) % 4;
               if ( ++count2 == 100 ) goto stuck;
               if ( base2 != base1 )
               {    if ( i == 0 || as_base(base2) != x[1][i-1] ) break;    }    }
          while(1)
          {    base3 = randomx( ) % 4;
               if ( ++count3 == 100 ) goto stuck;
               if ( base3 != base1 && base3 != base2 )
               {    if ( i == 0 || as_base(base3) != x[2][i-1] ) break;    }    }
          while(1)
          {    base4 = randomx( ) % 4;
               if ( ++count4 == 100 ) goto stuck;
               if ( base4 != base1 && base4 != base2 && base4 != base3 )
               {    if ( i == 0 || as_base(base4) != x[3][i-1] ) break;    }    }
          x[0] += as_base(base1);
          x[1] += as_base(base2);
          x[2] += as_base(base3);
          x[3] += as_base(base4);    }
     for ( int i = 0; i < 4; i++ )
          cout << x[i] << "\n";    }
