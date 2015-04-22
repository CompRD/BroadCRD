// Copyright (c) 2003 Broad Institute/Massachusetts Institute of Technology

#include "String.h"
#include "Vec.h"
#include "graphics/Color.h"
#include "graphics/HistogramsForDummies.h"
#include "graphics/PlotNormalHisto.h"
#include "random/NormalDistribution.h"



/*
 * RawPlot
 */
void RawPlot( const String &eps_file,
	      const String &title,
	      const vec<int> data,
	      float mean,
	      float stdev,
	      int max_bars,
	      bool plot_implied )
{
  using namespace ns_psplot;
  
  // Number of data points.
  longlong tot_count = 0;
  for (int ii=0; ii<(int)data.size( ); ii++)
    tot_count += data[ii];
  
  // Generate labels for the graph.
  String n_points = "Total data points = " + ToString( tot_count );
  String mean_stdev 
    = "mean = " + ToString( mean, 2 )
    + ", stdev = " + ToString( stdev, 2);
  
  vec<freetext> labels;
  labels.push_back( freetext( title, black, 16 ) );
  labels.push_back( freetext( n_points, black, 12 ) );
  labels.push_back( freetext( mean_stdev, black, 12 ) );

  // Parameters for BarPlot.
  color col = gray;
  color col_implied = red;
  float max_bar_height = 300;
  float bar_width = 6;
  float bar_sep = 1.4;
  short font_size = 7;
  float low = 0;
  int n_bars = Min( (int)data.size( ), max_bars );
  float high = n_bars - 1;
  vec<float> bot_tics;
  
  vec<float> f_data( n_bars );
  vec<float> f_implied( n_bars );
  for (int ii=0; ii<n_bars; ii++) {
    f_data[ii] = SafeQuotient( data[ii], tot_count );
    f_implied[ii] = NormalDensity( float( ii ), mean, stdev );
  }
  
  // Send graph to stream.
  ofstream f_fig( eps_file.c_str() );
  if ( plot_implied ) {
    bargraph graph = BiBarPlot( f_data, f_implied, low, high, bar_width,
				bar_sep, col, col_implied, font_size,
				max_bar_height, labels, true);
    f_fig << graph << "\n";
  }
  else {
    bargraph graph = BarPlot( f_data, low, high, bar_width, bar_sep,
			      col, font_size, max_bar_height, labels,
			      bot_tics, true);
    f_fig << graph << "\n";
  }
  f_fig.close( );
  
}



/*
 * BasicNormalHistoPlot
 */
void BasicNormalHistoPlot( const String &eps_file,
			   const String &title,
			   const vec<int> data,
			   int max_bars,
			   bool plot_implied )
{
  // Mean and standard deviation.
  longlong tot_count = 0;
  for (int ii=0; ii<(int)data.size( ); ii++)
    tot_count += data[ii];

  vec<float> data_points;
  data_points.reserve( tot_count );
  for (int ii=0; ii<(int)data.size( ); ii++)
    for (int jj=0; jj<data[ii]; jj++)
      data_points.push_back( ii );
  
  NormalDistribution nd = SafeMeanStdev( data_points );
  float mu = nd.mu_;
  float sigma = nd.sigma_;

  // Generate graph.
  RawPlot( eps_file, title, data, mu, sigma, max_bars, plot_implied );
}



/*
 * FatTailHistoPlot
 */
bool FatTailHistoPlot( const String &eps_file,
		       const String &title,
		       const vec<int> data,
		       int max_bars,
		       bool plot_implied )
{
  // Not enough data points.
  if ( data.size( ) < 8 )
    return false;
  
  // Locate the symmetry point in data. Skip the first entries.
  int symmetry = -1;
  for (int ii=3; ii<(int)data.size( )-1; ii++) {
    if ( data[ii] > data[ii+1] ) {
      symmetry = ii;
      break;
    }
  }
  if ( symmetry < 0 )
    return false;
  
  // Define a new set of data (mirror-image data on the right of symmetry).
  vec<int> sym_data( 2*symmetry + 1, 0 );
  sym_data[0] = 0;
  sym_data[symmetry] = data[symmetry];
  sym_data[sym_data.size( )-1] = 0;
  for (int ii=1; ii<symmetry; ii++) {
    sym_data[ii] = data[ii];
    sym_data[sym_data.size( ) - 1 - ii] = data[ii];
  }
  
  // Mean and standard deviation.
  longlong tot_count = 0;
  for (int ii=0; ii<(int)sym_data.size( ); ii++)
    tot_count += sym_data[ii];

  vec<float> data_points;
  data_points.reserve( tot_count );
  for (int ii=1; ii<(int)sym_data.size( )-1; ii++)
    for (int jj=0; jj<sym_data[ii]; jj++)
      data_points.push_back( ii );
  
  NormalDistribution nd = SafeMeanStdev( data_points );
  float mu = nd.mu_;
  float sigma = nd.sigma_;
  
  // Generate graph and return.
  RawPlot( eps_file, title, data, mu, sigma, max_bars, plot_implied );
  return true;
}



