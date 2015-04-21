/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "TokenizeString.h"
#include "math/Functions.h"
#include "random/NormalDistribution.h"
#include "reporting/SnapToGrid.h"
#include "util/RunCommand.h"

/**
 * PlotNormalLognormalFit
 *
 * Plot histogram and normal/log-normal Maximum Likelihood Estimates
 * of a list of numbers, loaded from a file (as floats). You need to
 * specify the range, and a size for the histograms' bins. The input
 * file should be a tab or space separated list of numbers, and you
 * can select which column to use.
 * 
 * IN: input file (as a .sample_seps.in from SamplePairedReadStats)
 * RANGE_LOW: low bin
 * RANGE_HIGH: high bin
 * DELTA: size of bins
 * COL: select colum, if the input file contains several columns (0-based)
 * GREP: if given, grep input lines
 * TMP: where tmp file will be saved (they are deleted on exit)
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( IN );
  CommandArgument_Int( RANGE_LOW );
  CommandArgument_Int( RANGE_HIGH );
  CommandArgument_Int_OrDefault( DELTA, 100 );
  CommandArgument_Int_OrDefault( COL, 0 );
  CommandArgument_String_OrDefault( GREP, "" );
  CommandArgument_String_OrDefault( TMP, "/tmp" );
  EndCommandArguments;

  // File names.
  String dataH_file = temp_file( TMP + "/MLE.dataH_XXXXXX" );
  String dataN_file = temp_file( TMP + "/MLE.dataN_XXXXXX" );
  String dataLN_file = temp_file( TMP + "/MLE.dataLN_XXXXXX" );
  String plot_file = temp_file( TMP + "/MLE.plot_XXXXXX" );

  // Load data, and sort it.
  vec<float> data;
  ifstream in( IN.c_str( ) );
  while ( 1 ) {
    String line;
    vec<String> tokens;
    getline( in, line );
    if ( ! in ) break;
    if ( GREP != "" && !line.Contains( GREP ) ) continue;
    Tokenize( line, tokens );
    int datum = tokens[COL].Int( );
    if ( datum < RANGE_LOW || datum > RANGE_HIGH ) continue;
    data.push_back( datum );
  }
  sort( data.begin( ), data.end( ) );
  
  // No data.
  if ( data.size( ) < 1 ) {
    cout << "No data found. Exit.\n" << endl;
    return 0;
  }

  // Generate log-data.
  float shift = 0;
  if ( data[0] < .0 ) shift = - data[0] + 1.0;
  for (int ii=0; ii<data.isize( ); ii++)
    data[ii] += shift;
  
  vec<float> logdata;
  logdata.reserve( data.size( ) );
  for (size_t ii=0; ii<data.size( ); ii++)
    logdata.push_back( log( data[ii] ) );
  
  // Compute normal and log-normal fits.
  NormalDistribution nd = SafeMeanStdev( data );
  NormalDistribution lognd = SafeMeanStdev( logdata );
  
  // The x axis, and the y axis.
  vec<float> x;
  int ii = 0;
  while ( 1 ) {
    float currx = RANGE_LOW + float( ( ii++ ) * DELTA );
    if ( currx > (float)RANGE_HIGH ) break;
    x.push_back( currx );
  }
  
  // The y axis.
  vec<float> yH;
  vec<float> yN;
  vec<float> yLN;

  {   // Histogram.
    yH.resize( x.size( ), 0 );
    for (int ii=0; ii<data.isize( ); ii++) {
      int bin = PlaceInBin( data[ii], x );
      if ( bin < 0 ) continue;
      yH[bin] += 1;
    }

    float bigsum = 0;
    for (int ii=0; ii<yH.isize( ); ii++)
      bigsum += yH[ii];

    ofstream out( dataH_file.c_str( ) );
    for (int ii=0; ii<yH.isize( ); ii++)
      out << x[ii] - shift << "\t" << yH[ii] / bigsum << "\n";
    out.close( );
  }

  {   // Plot of normal MLE.
    yN.resize( x.size( ), 0 );
    for (int ii=0; ii<x.isize( ); ii++)
      yN[ii] = DELTA * NormalDensity( x[ii], nd.mu_, nd.sigma_ );
    
    ofstream out( dataN_file.c_str( ) );
    for (int ii=0; ii<yN.isize( ); ii++)
      out << x[ii] - shift << "\t" << yN[ii] << "\n";
    out.close( );
  }

  {   // Plot of log-normal MLE.
    yLN.resize( x.size( ), 0 );
    for (int ii=0; ii<x.isize( ); ii++) {
      double amt = ( log( x[ii] ) - lognd.mu_ ) / ( lognd.sigma_ );
      yLN[ii]
	= DELTA
	* ( 1. / ( x[ii] * sqrt( 2.0 * M_PI ) * lognd.sigma_ ) ) 
	* exp( - 0.5 * amt * amt ) ;
    }

    ofstream out( dataLN_file.c_str( ) );
    for (int ii=0; ii<yLN.isize( ); ii++)
      out << x[ii] - shift << "\t" << yLN[ii] << "\n";
    out.close( );
  }

  // Gnuplot.
  ofstream oplot( plot_file.c_str( ) );

  String strNtitle = "normal fit (mean=" + ToString( nd.mu_ - shift, 1 ) + ")";
  
  String ln_median = ToString( exp( lognd.mu_ ) - shift, 1 );
  String strLNtitle = "log-normal fit (median=" + ln_median + ")";

  String str_sample = "(sample size: " + ToString( data.size( ) ) + " points)";
  String strH = + "using 1:2 title \"histogram\" with line,";
  String strN = "using 1:2 title \"" + strNtitle + "\" with line,";
  String strLN = "using 1:2 title \"" + strLNtitle + "\" with line";
  
  oplot << "set title \"" << Basename( IN ) << " " << str_sample << "\"\n"
	<< "\n"
	<< "plot \"" << dataH_file << "\" " << strH << " \\\n"
	<< "  \"" << dataN_file << "\" " << strN << " \\\n"
	<< "  \"" << dataLN_file << "\" " <<  strLN << "\n"
	<< endl;

  oplot.close( );
  
  // Run gnuplot.
  RunCommand( "gnuplot -persist " + plot_file );
  
  // Done.
  return 0;

}

