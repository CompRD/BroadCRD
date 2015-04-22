/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "random/NormalRandom.h"
#include "random/Poisson.h"

int PoissonKnuthCached( const double lambda, const vec<double>& unif, int& uptr )
{    long double L = expl(-lambda);
     if ( L > 0 )
     {    long double p = 1.0;
          int k = 0;
          do
          {    k++;
               p *= unif[uptr++];
               if ( uptr == unif.isize( ) ) uptr = 0;    }
          while ( p >= L );
          return k - 1;    }
     else return int( round( lambda + sqrt(lambda) * FastNormal( ) ) );    }

void PoissonRatioConfidenceInterval( const double lambda1, const double lambda2,
     const double p, const int sample, const vec<double>& unif, int& uptr, 
     double& low, double& high )
{    vec<double> x(sample);
     for ( int i = 0; i < sample; i++ )
     {    x[i] = double( PoissonKnuthCached( lambda1, unif, uptr ) )
               / double( PoissonKnuthCached( lambda2, unif, uptr ) );    }
     Sort(x);
     int n1 = int( round( (1.0-p)/2.0 * double(sample) ) );
     int n2 = int( round( ( 1.0 - (1.0-p)/2.0 ) * double(sample) ) );
     low = x[n1], high = x[n2];    }


