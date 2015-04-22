/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef POISSON_H
#define POISSON_H

#include "CoreTools.h"

// PoissonKnuthCached.  Generate Poisson-distributed random numbers.  This is a
// simple algorithm found on Wikipedia, apparently due to Knuth, adopted for speed
// to work off a fixed vector of uniformly distributed random numbers unif
// (0 <= unif[j] < 1) that could be computed e.g. with drand48( ).  The integer
// uptr points to an entry of unif.  It is successively advanced.
//
// Modification: If lambda is large enough to make expl(-lambda) = 0, which on one 
// architecture is ~11400, then we instead treat the distribution as normal 
// and compute using FastNormal.
//
// This algorithm is almost certainly not quite correct and not as fast as it
// could be.  There are several articles in the literature that provide a more
// sophisticated approach.

int PoissonKnuthCached( const double lambda, const vec<double>& unif, int& uptr );

// PoissonRatioConfidenceInterval.  Return a confidence interval [low,high] for 
// probability p (0 < p < 1) for the ratio of two Poisson-distributed random 
// variables with respective means lambda1 and lambda2.  The variables unif and 
// uptr are as in PoissonKnuthCached.
//
// I have no idea how this should really be computed.

void PoissonRatioConfidenceInterval( const double lambda1, const double lambda2, 
     const double p, const int sample, const vec<double>& unif, int& uptr,
     double& low, double& high );

#endif
