// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
//


/*
 * Cambridge, September 14, 2001
 *
 * HistogramsForDummies.cc
 *
 */
#include <math.h>

#include "graphics/BarGraph.h"
#include "graphics/BuildHistogram.h"
#include "math/Functions.h"
#include "graphics/HistogramsForDummies.h"
#include "String.h"
#include "system/System.h"
#include "Vec.h"



/*
 * Constants
 *
 * MIN_BARS:
 *  Minimum number of strictly > 0 bars in histogram.
 *
 * MAX_BARS:
 *  Maximum number of bars in histogram.
 *
 * MIN_BAR_HEIGHT:
 *  Minimum bar height: if none of the bar is at least this tall, then
 *  the graph is not produced.
 */
const int MIN_BARS = 5;
const int MAX_BARS = 55;
const int MIN_BAR_HEIGHT = 7;



/*
 * GenerateHistogram
 */
int GenerateHistogram( vec<ns_psplot::freetext> labels,
		       const vec<int>& lengths,
		       String file_name,
		       String file_data,
		       float cutoff,
		       float mean,
		       float stdev, 
                       const int digits_to_right )
{
  vec<float> f_lengths;
  f_lengths.reserve( lengths.size() );

  for (int ii=0; ii<(int)lengths.size(); ii++ )
    f_lengths.push_back( (float)lengths[ii] );

  return GenerateHistogram(labels, f_lengths, file_name, file_data, cutoff, mean, stdev, digits_to_right);
}



/*
 * GenerateHistogram
 */
int GenerateHistogram( vec<ns_psplot::freetext> labels,
		       const vec<float>& lengths,
		       String file_name,
		       String file_data,
		       float cutoff,
		       float mean,
		       float stdev, 
                       const int digits_to_right )
{
  using namespace ns_psplot;
  
  vec<float> loc_lengths;
  loc_lengths.reserve( lengths.size() );

  // If cutoff > 0, must eliminate any length which is more than
  // ( cutoff * stdev ) away from mean.
  if ( cutoff > 0 ) {
    // If the need be, calculates mean and stdev (can be expensive!)
    if ( 0 == mean || 0 == stdev ) {
      NormalDistribution moments = SafeMeanStdev( lengths );
      
      mean = moments.mu_;
      stdev = moments.sigma_;
    }

    int lev = 0;
    while ( lev < (int)lengths.size() ) {
      if ( fabs( lengths[lev] - mean ) < cutoff * stdev )
	loc_lengths.push_back( lengths[lev] );
      lev++;
    }
  }
  else {
    loc_lengths = lengths;
  }

  // Build actual histograms.
  vec<float> x_bars;
  vec<float> y_bars;

  // Build the histogram.
  int num_bars = 30;
  float_histogram histo( loc_lengths );
  histo.SetBarNumber( num_bars );
  int check_histo = histo.ProduceHistogram(x_bars, y_bars);

  // Error check.
  if ( 0 == check_histo )
    return 0;

  // We want at least MIN_BARS bars, and no more than MAX_BARS bars..
  if ( (int)y_bars.size() < MIN_BARS || (int)y_bars.size() > MAX_BARS )
    return 0;

  // If none of the bars is at least MIN_BAR_HEIGHT tall, or most
  // of the bars are 0, then return.
  int n_positive_bars = 0;
  bool is_tall_enough = false;
  for (int ii=0; ii<(int)y_bars.size(); ii++) {
    if ( y_bars[ii] > 0 )
      ++n_positive_bars;

    if ( y_bars[ii] >= MIN_BAR_HEIGHT )
      is_tall_enough = true;

    if ( is_tall_enough && n_positive_bars >= MIN_BARS )
      break;
  }
  
  if ( !is_tall_enough || n_positive_bars < MIN_BARS )
    return 0;

  // Chop out short bars at the beginning and at the end (this
  //  may reduce the number of bars.) Leave at most two short
  //  bars at either side.
  float tallest_bar = -1;
  for (int ii=0; ii<(int)y_bars.size(); ii++)
    tallest_bar = ( y_bars[ii] > tallest_bar ) ? y_bars[ii] : tallest_bar;
  float short_bar = tallest_bar / 100.0;

  while ( num_bars > 3 &&
	  y_bars[0] < short_bar &&
	  y_bars[1] < short_bar && 
	  y_bars[2] < short_bar ) {
    x_bars.erase( x_bars.begin() );
    y_bars.erase( y_bars.begin() );
    --num_bars;
  }

  while ( num_bars > 3 &&
	  y_bars[num_bars-1] < short_bar &&
	  y_bars[num_bars-2] < short_bar && 
	  y_bars[num_bars-3] < short_bar ) {
    x_bars.pop_back();
    y_bars.pop_back();
    --num_bars;
  }

  // Generate figure.
  //  max_bar_height: hight (in pixels) of the tallest bar;
  //  bar_width: width of each bar (in pixels;)
  //  bar_sep: distance between bars (in pixels;)
  //  col: color for the bars (see BarGraph.h;)
  //  font_size: font size for the tags (tic labels.)
  float max_bar_height = 150;
  float bar_width = 6;
  float bar_sep = 1.5;
  color col = gray;
  short font_size = 7;

  // Create the file.
  using namespace ns_psplot;

  file_name += ".eps";
  ofstream f_fig( file_name.c_str() );

  float low = x_bars[0];
  float high = x_bars[x_bars.size()-1];
  vec<float> bot_tics;

  bargraph graph = BarPlot(y_bars, low, high, bar_width, bar_sep,
			   col, font_size, max_bar_height, labels,
			   bot_tics, true, digits_to_right );
  
  f_fig << graph << "\n";

  // Create data file.
  if ( "" != file_data ) {
    ofstream out( file_data.c_str() );
    for (int ii=0; ii<(int)x_bars.size(); ii++)
      out << x_bars[ii] << "\t" << y_bars[ii] << "\n";
  }
  
  //  All right, can return.
  return 1;
}



