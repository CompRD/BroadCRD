/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include <math.h>

#include "system/Assert.h"
#include "math/Functions.h"
#include "math/LinMinSolve.h"
#include "math/Matrix.h"
#include "random/Random.h"
#include "system/System.h"
#include "system/Types.h"
#include "Vec.h"

// MaxAbs: return the maximum of the absolute values of the entries of a vector.

template<class T> T MaxAbs( const vec<T>& v )
{    ForceAssert( v.size( ) > 0 );
     T answer = 0;
     for ( int i = 0; i < (int) v.size( ); i++ )
          answer = Max( answer, Abs( v[i] ) );
     return answer;    }

// DiscrepToSolveNonneg: Given a matrix A and a vector b, try to approximate the 
// nonnegative vector x which minimizes max|Ax-b|.  Return max|Ax-b| for this x.
//
// (Note: In the above notation, if v is a vector, |v| is a vector of the same
// size, with each of its entries replaced by its absolute value.  If w is a
// vector, max w is the maximum of its entries.)
//
// Method: iterative.  Make a small random perturbation p of x.  If it is an
// improvement, find a scalar m which optimizes x + mp.  Replace x by it.  Then
// go back and look for another p.
//
// This method is -- to a certain extent -- justified by the following observation:
// if x minimizes max|Ax-b|, and y >= 0, then along the path from y to x, the 
// value of w |--> max|Aw-b| is decreasing.  Therefore at any point which is not
// optimal, any neighborhood of it contains a better point.
//
// The code contains heuristic constants which are not universal.  They just
// happen to work for the types of problems we've run this on.
//
// There must be better methods.  This method is slow and sometimes (rarely) gives 
// the wrong answer.
//
// We give some special (exact) solutions because they are fast and needed for
// the applications.

float DiscrepToSolveNonneg( const matrix<float>& A, const vec<float>& b,
     Bool no_special )
{    int r = A.Nrows( ), c = A.Ncols( );
     ForceAssert( r == (int) b.size( ) );
     ForceAssert( r > 0 );
     ForceAssert( c > 0 );
               
     // Give a special solution in the case where r = 1 and c = 1.

     if ( r == 1 && c == 1 )
     {    if ( A[0][0] == 0 ) return Abs( b[0] );
          float x = b[0] / A[0][0];
          if ( x >= 0 ) return 0.0;
          return Abs( b[0] );    }

     // Give a special solution in the case where r = 2 and c = 1.

     if ( r == 2 && c == 1 )
     {    
          float a1 = A[0][0], a2 = A[1][0];
          float b1 = b[0], b2 = b[1];
          if ( a1 < 0 )
          {    a1 = -a1;
               b1 = -b1;    }
          if ( a2 < 0 )
          {    a2 = -a2;
               b2 = -b2;    }

          // Looking at a1x - b1, a2x - b2, where a1 >= 0 and a2 >= 0.

          if ( a1 == 0 && a2 == 0 ) return MaxAbs(b);
          if ( a1 == 0 )
          {    float x = b2/a2;
               if ( x <= 0 ) return MaxAbs(b);
               else return Abs(b1);    }
          if ( a2 == 0 )
          {    float x = b1/a1;
               if ( x <= 0 ) return MaxAbs(b);
               else return Abs(b2);    }
          float x1 = b1/a1, x2 = b2/a2;
          if ( x1 <= 0 && x2 <= 0 ) return MaxAbs(b);
          else if ( x1 > 0 && x2 <= 0 )
          {    float x = (b1+b2)/(a1+a2);
               if ( x > 0 ) return Abs( a1*x - b1 );
               else return MaxAbs(b);    }
          else if ( x1 <= 0 && x2 > 0 )
          {    float x = (b1+b2)/(a1+a2);
               if ( x > 0 ) return Abs( a2*x - b2 );
               else return MaxAbs(b);    }
          else return Abs( a1*(b1+b2)/(a1+a2) - b1 );    }

     // Give a special solution in the case where c = 1.

     if ( c == 1 )
     {    static vec<float> candidates;
          candidates.clear( );
          candidates.push_back(0);
          for ( int j = 0; j < r; j++ )
               if ( A[j][0] != 0 ) candidates.push_back( b[j] / A[j][0] );
          for ( int i = 0; i < r; i++ )
               for ( int j = i + 1; j < r; j++ )
               {    if ( A[i][0] - A[j][0] != 0 )
                         candidates.push_back( (b[i]-b[j]) / (A[i][0] - A[j][0]) );
                    if ( -A[i][0] - A[j][0] != 0 )
                         candidates.push_back( (-b[i]-b[j]) / (-A[i][0] - A[j][0]) );
                    if ( A[i][0] + A[j][0] != 0 )
                         candidates.push_back( (b[i]+b[j]) / (A[i][0] + A[j][0]) );
                    if ( -A[i][0] + A[j][0] != 0 )
                         candidates.push_back( (-b[i]+b[j]) / (-A[i][0] + A[j][0]) );
                         }
          float best_val = 0;
          for ( int i = 0; i < (int) candidates.size( ); i++ )
          {    float x = candidates[i];
               if ( x >= 0 )
               {    float val = 0;
                    for ( int j = 0; j < r; j++ )
                         val = Max( val, Abs( A[j][0] * x - b[j] ) );
                    if ( i == 0 ) best_val = val;
                    else if ( val < best_val ) best_val = val;    }    }
          return best_val;    }

     // Give the solution in a very special case.

     if ( !no_special && r == 2 && c == 2 && A[0][0] > 0 && A[0][1] == 0
          && A[1][0] > 0 && A[1][0] == A[1][1] )
     {
          // cout << "in very special case" << endl; // XXX
          // Case 1: optimum is where x[0] = 0 and x[1] = 0.

          float case1 = MaxAbs(b);

          // Case 2: optimum is where x[0] > 0 and x[1] = 0.

          static matrix<float> As(2, 1);
          As[0][0] = A[0][0];
          As[1][0] = A[1][0];
          float case2 = DiscrepToSolveNonneg( As, b );

          // Case 3: optimum is where x[0] = 0 and x[1] > 0.

          As[0][0] = A[0][1];
          As[1][0] = A[1][1];
          float case3 = DiscrepToSolveNonneg( As, b );

          // Case 4: optimum is where x[0] > 0 and x[1] > 0.  This is 
          // only possible in the case of a perfect solution.  (Up to
          // this case, the only hypotheses we used were that r = 2 and c = 2.)

          float x1 = b[0] / A[0][0];
          float x2 = (b[1] - A[1][0] * x1) / A[1][1];

          if ( x1 >= 0 && x2 >= 0 ) return 0;
          else return std::min( {case1, case2, case3} );    }

     // Now do the general case.

     static vec<float> x, x_alt, diff, p;
     x.resize_and_set( c, 0 );
     x_alt.resize(c);
     p.resize(c);
     diff.resize(r);
     mulsub( A, x, b, diff );
     float best = MaxAbs(diff);

     // Heuristic constants, which should really depend on the sort of data
     // given as input and the goal.  There are other heuristic constants,
     // embedded in the code.

     const float eps = 0.000001;
     const int delta = 10;
     const int tries = 500;
     const int max_improve_dirs = 40;

     int improve_dirs = 0;
     for ( int k = 0; k < tries; k++ )
     {    
          for ( int j = 0; j < c; j++ )
          {    p[j] = (randomx( ) % (2 * delta)) - delta;
               if ( x[j] < eps && p[j] < 0 ) p[j] = 0;    }

          for ( int j = 0; j < c; j++ )
          {    x_alt[j] = x[j] + p[j];
               if ( x_alt[j] < eps ) x_alt[j] = 0;    }
          mulsub( A, x_alt, b, diff );
          float maybe_best = MaxAbs(diff);
          if ( maybe_best < best )
          {    ++improve_dirs;
               static vec<float> q, a;
               q.resize(r);
               a.resize(r);
               mulsub( A, x, b, q );
               for ( int i = 0; i < r; i++ )
                    q[i] = -q[i];
               mul( A, p, a );
               static vec<float> candidates;
               candidates.clear( );
               candidates.push_back(0);
               for ( int j = 0; j < r; j++ )
                    if ( a[j] != 0 ) candidates.push_back( q[j] / a[j] );
               for ( int i = 0; i < r; i++ )
                    for ( int j = i + 1; j < r; j++ )
                    {    if ( a[i] - a[j] != 0 )
                              candidates.push_back( (q[i]-q[j]) / (a[i] - a[j]) );
                         if ( -a[i] - a[j] != 0 )
                              candidates.push_back( (-q[i]-q[j]) / (-a[i] - a[j]) );
                         if ( a[i] + a[j] != 0 )
                              candidates.push_back( (q[i]+q[j]) / (a[i] + a[j]) );
                         if ( -a[i] + a[j] != 0 )
                              candidates.push_back( (-q[i]+q[j]) / (-a[i] + a[j]) );
                              }
               for ( int i = 0; i < c; i++ )
                    if ( p[i] != 0 ) candidates.push_back( -x[i] / p[i] );
               float best_val = 0;
               int best_candidate = 0;
               for ( int i = 0; i < (int) candidates.size( ); i++ )
               {    float t = candidates[i];
                    if ( t >= 0 )
                    {    int j;
                         for ( j = 0; j < c; j++ )
                         {    x_alt[j] = x[j] + t * p[j];
                              if ( x_alt[j] < 0 ) break;    }
                         if ( j < c ) continue;
                         float val = 0;
                         for ( j = 0; j < r; j++ )
                              val = Max( val, Abs( a[j] * t - q[j] ) );
                         if ( i == 0 ) best_val = val;
                         else if ( val < best_val ) 
                         {    best_val = val;    
                              best_candidate = i;    }    }    }
               if ( best_val / best < 0.99999 ) // was 9999
               {    improve_dirs = 0;
                    k = 0;    }
               best = best_val;
               // PRINT(best); // XXX
               float t = candidates[best_candidate];
               for ( int j = 0; j < c; j++ )
                    x[j] += t * p[j];
               if ( improve_dirs == max_improve_dirs ) break;    }    }
               
     return best;    }
