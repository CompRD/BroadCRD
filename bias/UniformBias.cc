/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "bias/UniformBias.h"
#include "math/Functions.h"
#include "random/Random.h"

void Smooth( const vec< vec<int> >& x, int radius, vec< vec<int> >& s )
{    ForceAssertGe( radius, 1 );
     --radius;
     s.resize( x.size( ) );
     double dr = double( 2*radius + 1 );
     for ( int c = 0; c < x.isize( ); c++ )
     {    s[c].resize( x[c].size( ) );
          ForceAssert( x[c].nonempty( ) );
          longlong sum = x[c][0], denom = 1;
          for ( int j = 1; j <= radius; j++ )
          {    if ( j < x[c].isize( ) )
               {    sum += x[c][j];
                    ++denom;    }    }
          s[c][0] = int( round( ( dr * double(sum) ) / double(denom) ) );
          longlong last_sum = sum, last_denom = denom;
          for ( int i = 1; i < x[c].isize( ); i++ )
          {    longlong sum = last_sum, denom = last_denom;
               if ( i - radius > 0 ) sum -= x[c][i-radius-1];
               else ++denom;
               if ( i + radius < x[c].isize( ) ) sum += x[c][i+radius];
               else --denom;
               if ( denom == 2*radius + 1 ) s[c][i] = sum;
               else
               {    s[c][i] = int( round( ( dr * double(sum) )
                         / double(denom) ) );    }    
               last_sum = sum, last_denom = denom;    }    }    }

double HowBiased( const vec< vec<int> >& x, const vec<int>& N )
{    double B = 0.0;
     longlong ss = 0;
     for ( int c = 0; c < x.isize( ); c++ )
     {    for ( int i = 0; i < x[c].isize( ); i++ )
               ss += x[c][i];    }
     double u = double(ss) / double( BigSum(N) );
     for ( int c = 0; c < x.isize( ); c++ )
     {    double sum = 0.0;
          for ( int i = 0; i < N[c]; i++ )
          {    double y = double( x[c][i] ) - u;
               sum += y * y;    }
          B += double(N[c]) 
               * 100.0 * ( sum / ( u * u * double(N[c]) ) - 1.0/u );    }
     return B / double( BigSum(N) );    }

double HowBiasedRelRaw( const vec< vec<int> >& x1, const vec< vec<int> >& x2,
     const vec<int>& N )
{    ForceAssertEq( x1.size( ), x2.size( ) );
     longlong ss = 0, ssr = 0, x1sum = 0, x2sum = 0, nc = 0;
     for ( int c = 0; c < x1.isize( ); c++ )
     {    ForceAssertEq( x1[c].size( ), x2[c].size( ) );
          nc += x1[c].size( );
          for ( int i = 0; i < x1[c].isize( ); i++ )
          {    x1sum += x1[c][i], x2sum += x2[c][i];
               ss += x1[c][i] * x2[c][i];    }    }
     vec<int> x1r( nc, 0 ), x2r( nc, 0 );
     for ( int i = 0; i < x1sum; i++ )
          ++x1r[ randomx( ) % nc ];
     for ( int i = 0; i < x2sum; i++ )
          ++x2r[ randomx( ) % nc ];
     for ( int i = 0; i < nc; i++ )
          ssr += x1r[i] * x2r[i];
     // return 20000.0 * double(ss-ssr) / double(x1sum + x2sum);    }
     return 4.0 * 100.0 
          * double(BigSum(N)) * double(ss-ssr) / pow(x1sum + x2sum, 2.0);    }

double HowBiasedDev( const vec< vec<int> >& x, const vec<int>& N )
{    double D = 0.0;
     longlong ss = 0;
     for ( int c = 0; c < x.isize( ); c++ )
     {    for ( int i = 0; i < x[c].isize( ); i++ )
               ss += x[c][i];    }
     double u = double(ss) / double( BigSum(N) );
     for ( int c = 0; c < x.isize( ); c++ )
          D += double(N[c]) * u * u;
     return 100.0 * sqrt( 2.0 / D );    }
