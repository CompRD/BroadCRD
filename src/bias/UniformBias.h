/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// HowBiased[Dev].  Given points on one or more intervals, return a metric 
// (HowBiased) that assesses the extent to which the points land uniformly on the 
// intervals, and error bars (+/- HowBiasedDev) for that metric.  The metric is 
// INTENDED to have the following APPROXIMATE properties:
// - It is additive, i.e. it returns zero for randomly distributed points, and for
// two successive processes that introduce small amounts of bias, the value of the 
// metric for the output of the joint process is approximately the sum of the values 
// for the outputs of the two processes taken separately.
// - It is independent of the point density.
// - It is independent of how long the intervals are.
// - It will return a value of 100 if half of the interval space has zero coverage,
// and the other half has random coverage.
// - It will return a negative value if the points are distributed more uniformly
// than random.

#ifndef UNIFORM_BIAS_H
#define UNIFORM_BIAS_H

#include "CoreTools.h"

// HowBiased: given points x1 from [0,N), assign a bias score,
// 100 * [ (1/N) * integral_0^N ( c(x)/u - 1 )^2 dx - 1/u ],
// were u is the mean coverage and c(x) is the number of points at x.
// This is just an appropriately normalized version of the variance, designed to
// assess the lack of uniform coverage over and above what's expected if the data is
// truly random.
//
// This generalizes to the cases where one has points x1 from [0,N1), 
// x2 from [0,N2), etc.

double HowBiased( const vec< vec<int> >& x, const vec<int>& N );

// HowBiasedDev: in terms of the notation of HowBiased, 
// return 100 * sqrt( 2 / ( N * u^2 ) ).

double HowBiasedDev( const vec< vec<int> >& x, const vec<int>& N );

double HowBiasedRelRaw( const vec< vec<int> >& x1, const vec< vec<int> >& x2,
     const vec<int>& N );

// Smooth: replace x[c][i] by the sum of x[c][j], i-radius < j < i+radius.

void Smooth( const vec< vec<int> >& x, int radius, vec< vec<int> >& s );

#endif
