// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
//

#include "MainTools.h"
#include "system/ParsedArgs.h"
#include "graphics/BarGraph.h"
#include "graphics/BuildHistogram.h"
#include "lookup/LookAlign.h"
#include "lookup/PrintLookAlign.h"



/*
 * ErrorRatesHistogram
 *
 * Load the given look_aligns and plot the histogram of the error
 * rates (as reported in look_aligns). No filtering is performed. It
 * saves the histogram both as text file and as eps, in
 * OUT_DIR/OUT_BASE.txt and OUT_DIR/OUT_BASE.eps.
 *
 * HITS_FILE: full path name for look_aligns (input)
 * OUT_DIR: full path name for output dir
 * OUT_BASE: base name for output files
 * INCREMENT: for bins in histogram
 * MAX_BIN: largest bin to plot
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( HITS_FILE );
  CommandArgument_String( OUT_DIR );
  CommandArgument_String_OrDefault( OUT_BASE, "histogram" );
  CommandArgument_Double_OrDefault( INCREMENT, 0.01 );
  CommandArgument_Double_OrDefault( MAX_BIN, 1.1 );
  EndCommandArguments;
  
  // Dir and file names.
  String out_txt_file = OUT_DIR + "/" + OUT_BASE + ".txt";
  String out_eps_file = OUT_DIR + "/" + OUT_BASE + ".eps";
  
  Mkdir777( OUT_DIR );

  // Load.
  cout << Date( ) << ": loading hits" << endl;
  vec<look_align_plus> hits;
  LoadLookAlignPlus( HITS_FILE, hits );
  
  // Generate vector with error rates.
  cout << Date( ) << ": generate histogram" << endl;
  vec<float> rates;
  rates.reserve( hits.size( ) );
  for (int ii=0; ii<(int)hits.size( ); ii++)
    rates.push_back( hits[ii].ErrorRate( ) );
  sort( rates.begin( ), rates.end( ) );
  
  // Generate histogram.
  float increment = INCREMENT;
  vec<float> x_bars( 1, 0 );
  while ( x_bars.back( ) < float( MAX_BIN ) )
    x_bars.push_back( (float)x_bars.size( ) * increment );
  
  vec<float> y_bars( x_bars.size( ), 0 );
  
  float_histogram histo( rates );
  if ( 0 == histo.ProduceHistogram( x_bars, y_bars, false ) ) {
    cout << "\n Failed to generate histogram. Exit\n" << endl;
    return 0;
  }

  // Save txt.
  cout << Date( ) << ": saving in " << OUT_DIR << endl;
  ofstream out_txt( out_txt_file.c_str( ) );
  out_txt << "Error Rate Histogram for\n" << HITS_FILE << "\n";
  for (int ii=0; ii<(int)x_bars.size( ); ii++)
    out_txt << x_bars[ii] << "\t" << y_bars[ii] << "\n";
  out_txt.close( );

  // Generate and save eps.
  using namespace ns_psplot;

  float max_bar_height = 150;
  float bar_width = 6;
  float bar_sep = 1.5;
  color col_label = black;
  color col_bars = gray;
  short font_size = 7;

  vec<freetext> labels( 2 );
  labels[0] = freetext( "Error Rate Histogram for", col_label, 12 );
  labels[1] = freetext( HITS_FILE, col_label, 10 );

  float low = *min_element( x_bars.begin( ), x_bars.end( ) );
  float high = *max_element( x_bars.begin( ), x_bars.end( ) );
  vec<float> bot_tics;

  bargraph graph = BarPlot( y_bars, low, high, bar_width, bar_sep,
			    col_bars, font_size, max_bar_height, labels,
			    bot_tics, true);
  
  ofstream out_eps( out_eps_file.c_str() );
  out_eps << graph << "\n";
  out_eps.close( );
  
  // Done.
  cout << Date( ) << ": done" << endl;

}
