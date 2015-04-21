/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Quality.h"
#include "MainTools.h"
#include "String.h"
#include "graphics/BarGraph.h"

/**
 * QualsHistogram
 *
 * Load the given vecqualvector and generate histogram and distribution
 * for the quality scores (both as eps and text).
 *
 * QUALS: input vecqualvector, full path name
 * OUTDIR: full path name for output dir
 * CAP: cap quality scores
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( QUALS );
  CommandArgument_String( OUTDIR );
  CommandArgument_UnsignedInt_OrDefault( CAP, 50 );
  EndCommandArguments;

  // File names.
  String raw_file = OUTDIR + "/raw.data";
  String histo_eps_file = OUTDIR + "/histo.eps";
  String histo_txt_file = OUTDIR + "/histo.txt";
  String dist_eps_file = OUTDIR + "/distrib.eps";
  String dist_txt_file = OUTDIR + "/distrib.txt";
  Mkdir777( OUTDIR );

  // Load.
  cout << Date( ) << ": loading quals" << endl;
  vecqualvector quals;
  quals.ReadAll( QUALS );

  // Labels.
  String str_histo = "Histogram of quality scores for";
  String str_dist = "Distribution of quality scores for";

  // Generate raw data.
  cout << Date( ) << ": parsing " << quals.size( ) << " quals" << endl;
  vec<longlong> rawbins( 1 + (int)CAP, 0 );
  for (vecqvec::size_type ii=0; ii<quals.size( ); ii++)
    for (qvec::size_type jj=0; jj<quals[ii].size( ); jj++)
      rawbins[ Min( (int)quals[ii][jj], (int)CAP ) ] += 1;

  // Save raw data to file.
  longlong sum = 0;
  for (int ii=0; ii<(int)rawbins.size( ); ii++)
    sum += rawbins[ii];

  ofstream raw_out( raw_file.c_str( ) );
  raw_out << "Raw data for" << QUALS << "\n\n"
	  << "Total bases: " << sum << "\n\n";
  for (int ii=0; ii<(int)rawbins.size( ); ii++)
    raw_out << ii << "\t" << rawbins[ii] << "\n";
  raw_out.close( );

  // Normalize.
  cout << Date( ) << ": generating histogram" << endl;
  vec<float> bins( rawbins.size( ), 0 );
  for (int ii=0; ii<(int)bins.size( ); ii++)
    bins[ii] = 100.0 * SafeQuotient( rawbins[ii], sum );

  // Generate figure.
  using namespace ns_psplot;

  float low = 0;
  float high = Max( bins );
  vec<float> bot_tics( 1 + (int)CAP, 0 );
  for (int ii=0; ii<(int)bot_tics.size( ); ii++)
    bot_tics[ii] = (float)ii;

  float max_bar_height = 150;
  float bar_width = 6;
  float bar_sep = 1.5;
  short font_size = 7;
  color col = gray;

  vec<ns_psplot::freetext> h_labels( 2 );
  h_labels[0] = freetext( str_histo, black, 12 );
  h_labels[1] = freetext( QUALS, black, 9 );

  bargraph h_graph = BarPlot( bins, low, high, bar_width, bar_sep,
			      col, font_size, max_bar_height, h_labels,
			      bot_tics, true);

  // Save histogram.
  cout << Date( ) << ": saving histogram" << endl;

  ofstream histo_eps_out( histo_eps_file.c_str( ) );
  histo_eps_out << h_graph << "\n";
  histo_eps_out.close( );

  ofstream histo_txt_out( histo_txt_file.c_str( ) );
  histo_txt_out << str_histo << "\n" << QUALS << "\n\n";
  for (int ii=0; ii<(int)bins.size( ); ii++)
    histo_txt_out << ii << "\t" << ToString( bins[ii], 2 ) << "%\n";
  histo_txt_out.close( );

  // Generate distribution
  cout << Date( ) << ": generating distribution" << endl;
  for (int ii=(int)bins.size( )-2; ii>=0; ii--)
    bins[ii] = bins[ii] + bins[ii+1];

  // Generate figure.
  high = Max( bins );

  vec<ns_psplot::freetext> d_labels( 2 );
  d_labels[0] = freetext( str_dist, black, 12 );
  d_labels[1] = freetext( QUALS, black, 9 );

  bargraph d_graph = BarPlot( bins, low, high, bar_width, bar_sep,
			      col, font_size, max_bar_height, d_labels,
			      bot_tics, true);

  // Save distribution .
  cout << Date( ) << ": saving distribution" << endl;

  ofstream dist_eps_out( dist_eps_file.c_str( ) );
  dist_eps_out << d_graph << "\n";
  dist_eps_out.close( );

  ofstream dist_txt_out( dist_txt_file.c_str( ) );
  dist_txt_out << str_dist << "\n" << QUALS << "\n\n";
  for (int ii=0; ii<(int)bins.size( ); ii++)
    dist_txt_out << ii << "\t" << ToString( bins[ii], 2 ) << "%\n";
  dist_txt_out.close( );

  // Done.
  cout << Date( ) << ": done" << endl;
}
