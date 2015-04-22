/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "math/Functions.h"
#include "String.h"
#include "Vec.h"
#include "graphics/CContigPlotter.h"
#include "graphics/PlotCoverage.h"
#include "util/PickTic.h"

/**
 * PlotCoverage
 */
void PlotCoverage( int win_size,
		   const vec<float> &cov,
		   const String &outfile )
{
  float mult = 2.5;

  int length = cov.size( );
  int hop = win_size / 2;

  vec<float> all_averages;

  CDataPlotter data( red );
  int data_size = 0;
  while ( data_size * hop < (int)cov.size( ) )
    data_size++;
  data.Reserve( data_size );
  if ( mult > 0 )
    all_averages.reserve( data_size );
  int post = 0;
  float data_max = 0;
  while ( post * hop < (int)cov.size( ) ) {
    longlong totcov = 0;
    int begin = post * hop;
    int end = Min( begin + win_size, (int)cov.size( ) );
    float interval_size = end - begin;
    for (int ii=begin; ii<end; ii++)
      totcov += (int)cov[ii];
    float average = float( totcov ) / interval_size;
    if ( mult > 0 )
      all_averages.push_back( average );
    
    data.AddPoint( begin, average );
    data_max = Max( data_max, average );
    post++;
  }
  
  CDataPlotter median( blue );
  CDataPlotter median_high ( blue );
  if ( mult > 0 ) {
    sort( all_averages.begin( ), all_averages.end( ) );
    float med = Median( all_averages, 0 );
    float med_h = mult * med;
    median.AddPoint( 0, med );
    median.AddPoint( length, med );
    median_high.AddPoint( 0, med_h );
    median_high.AddPoint( length, med_h );
  }

  vec<CTic> xtics;
  {
    float max = length;
    float unit = 1000000;
    int power = 0;
    int max_tics = 30;
    PickTic( max, unit, max_tics, power );
    float step = pow( 10.0, power ) * unit;
    float pos = 0;
    while ( pos < length ) {
      String name = ToString( pos / unit , 2 ) + ( pos == 0 ? " Mb" : "" );
      CTic tic( pos, name );
      xtics.push_back( tic );
      pos += step;
    }
  }
  
  vec<CTic> ytics;
  {
    float max = data_max;
    float unit = 1.0;
    int power = 0;
    int max_tics = 12;
    if ( max > 0 )
      PickTic( max, unit, max_tics, power );
    float step = pow( 10.0, power ) * unit;
    float pos = 0;
    while ( pos < max ) {
      String name = ToString( pos / unit , 1 ) + "X";
      CTic tic( pos, name );
      ytics.push_back( tic );
      pos += step;
    }
  }
  
  CContigPlotter plotter( length );
  if ( mult > 0 ) {
    plotter.AddData( median );
    plotter.AddData( median_high );
  }
  plotter.AddData( data );
  plotter.AddHorizontalTics( xtics );
  plotter.AddVerticalTics( ytics );
  
  plotter.GoFigure( outfile );
}

/**
 * PlotCoverage
 */
bool PlotCoverage( int target_id,
		   int target_len,
		   int window_size,
		   const vec<seq_interval> &coverages,
		   const vec<int> &fcov,
		   const String &eps_file )
{
  // Test for early exit.
  int firstcov = fcov[ target_id ];
  if ( firstcov < 0 || target_len < window_size )
    return false;
  
  // The vector with base coverage for target_id.
  vec<float> raw( target_len, 0 );
  for (int cov_id=firstcov; cov_id<(int)coverages.size( ); cov_id++) {
    if ( coverages[cov_id].SeqId( ) != target_id )
      break;
    int begin = coverages[cov_id].Begin( );
    int end = coverages[cov_id].End( ); 
    int amount = coverages[cov_id].IntervalId( );
    for (int jj=begin; jj<end; jj++) {
      raw[jj] = (float)amount;
    }
  }
  
  // Plot this target coverage and return.
  PlotCoverage( window_size, raw, eps_file );
  return true;
}

/**
 * PlotCoverage
 */
bool PlotCoverage( int window_size,
		   const qualvector &cov,
		   const String &eps_file )
{
  // Test for early exit.
  if ( (int)cov.size( ) < window_size )
    return false;
  
  // The vector with base coverage for target_id.
  vec<float> raw( cov.size( ), 0 );
  for (int ii=0; ii<(int)cov.size( ); ii++)
    raw[ii] = (float)(int)cov[ii];
  
  // Plot this target coverage and return.
  PlotCoverage( window_size, raw, eps_file );
  return true;
}

