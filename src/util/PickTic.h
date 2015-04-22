/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef PICK_TIC_H
#define PICK_TIC_H

#include "math/Functions.h"
#include "String.h"
#include "graphics/CTic.h"
#include "graphics/Whiteboard.h"

/**
 * FitsTic
 */
inline bool FitsTic( const int max_tics,
		     const int power,
		     const float max,
		     const float unit )
{
  return ( unit * pow( 10.0, power ) * float( max_tics ) < max );
}

/**
 * PickTic
 *
 * We are looking at a plot of the real interval [0, max), and we want
 * to add tics to this segment in such a way that there are no more than
 * max_tics tics overall. Output is saved in power, and it is the "best"
 * choice for
 *
 *   tics_n =  n * ( unit * 10^power ); with n=0, 1, ...
 *
 * In other words, n is the largest possible integer (may be < 0) such
 * that there are less than max_tics tics of the type tics_n. Notice: 
 * max must be strictly positive.
 */
inline void PickTic( const float max,
		     const float unit,
		     const int max_tics,
		     int &power )
{
  ForceAssert( max > 0 );

  power = 0;
  while ( ! FitsTic( max_tics, power, max, unit ) )
    power--;
  while( FitsTic( max_tics, power, max, unit ) )
    power++;
}

/**
 * StockXTics
 *
 * Generate "stock" xtics (with unit=1000000). Result is stored in xtics.
 */
inline void StockXTics( vec<CTic> &xtics,
			const float max,
			const String &label = " Mb" )
{
  int power = 0;
  float unit = 1000000;
  int max_tics = 30;
  PickTic( max, unit, max_tics, power );
  float step = pow( 10.0, power ) * unit;
  float pos = 0;
  while ( pos < max ) {
    String name = ToString( pos / unit , 2 ) + ( pos == 0 ? label : "" );
    CTic tic( pos, name );
    xtics.push_back( tic );
    pos += step;
  }
}

/**
 * StockYTics
 *
 * Generate "stock" ytics (with unit=1). Result is stored in ytics.
 */
inline void StockYTics( vec<CTic> &ytics,
			const float max,
			const String &label = "" )
{
  int power = 0;
  float unit = 1.0;
  int max_tics = 12;
  if ( max > 0 )
    PickTic( max, unit, max_tics, power );
  float step = pow( 10.0, power ) * unit;
  float pos = 0;
  while ( pos < max ) {
    String name = ToString( pos / unit , 1 ) + ( pos == 0 ? label : "" );
    CTic tic( pos, name );
    ytics.push_back( tic );
    pos += step;
  }
}

#endif // PICK_TIC_H
