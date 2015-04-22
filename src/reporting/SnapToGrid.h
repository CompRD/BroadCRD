/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef SNAP_TO_GRID_H
#define SNAP_TO_GRID_H

#include "math/Functions.h"
#include "Vec.h"

/**
 * SnapToGrid
 *
 * Place a given double into the closest clustering bin. For example
 * if bins = { 0, 2.5, 4.5, 10, 40 } and x = 3.75 then x belongs to
 * bins[1], and SnapToGrid returns 1. If d(x,bins[i]) = d(x,bins[j])
 * then x can be placed either in bins[i] or bins[j]. On error it
 * returns -1. Which proves that reading the code is faster than
 * reading the documentation.
 */
inline int SnapToGrid( double x_val, const vec<double> &bins )
{
  if ( bins.size( ) < 1 )
    return -1;

  int pos = 0;
  double dist = Abs( bins[0] - x_val );
  for (int ii=0; ii<(int)bins.size( ); ii++) {
    double curr_dist = Abs( bins[ii] - x_val );
    if ( curr_dist < dist ) {
      pos = ii;
      dist = curr_dist;
    }
  }
  
  return pos;
}

/**
 * PlaceInBin
 *
 * Similar to SnapToGrid, but the snapping here is slightly different.
 * There are n bins, b(0), b(1), ..., b(n-1). With the convention that
 * b(n) = + \infty, PlaceInBin will place x_val in bin 0 if x_val < b(0);
 * or in bin ii if b(ii) <= x_val < b(ii+1), for ii={0, 1, ...,n-1}.
 * It assumes (but does not check) bins is non empty and sorted.
 */
template <class T>
inline int PlaceInBin( T x_val, const vec<T> &bins )
{
  typename vec<T>::const_iterator it;

  it = upper_bound( bins.begin( ), bins.end( ), x_val );
  if ( it == bins.begin( ) )
    return 0;
  return distance( bins.begin( ), it ) - 1;
}

#endif
