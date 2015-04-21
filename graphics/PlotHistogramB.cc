/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "FastIfstream.h"
#include "Histogram.h"
#include "TokenizeString.h"
#include "VecUtilities.h"

/**
 * PlotHistogramB
 *
 * It generates both text and an eps output files: <DATA>.{eps,histo}.
 *
 * DATA: can be a file or stdin (int only)
 * COL: print an histogram of the entries at the given column (0-based)
 */
int main( int argc, char *argv[] ) 
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( DATA );
  CommandArgument_Int_OrDefault( COL, 0 );
  
  // These are used to build histograms.
  const float MAX_MULT = 5.0;
  const vec<float> PROGRESSION = MkVec( float(25), float(50), float(100) );
  const int MIN_BINS = 60;
  const int MAX_BINS = 120;

  // File names.
  String eps_file = DATA + ".eps";
  String histo_file = DATA + ".histo";

  // Load data.
  vec<float> lens;

  istream *data;
  if ( DATA == "-" ) data = &cin ;
  else data = new ifstream( DATA.c_str( ) );

  String line;
  vec<String> tokens;
  while( 1 ) {
    getline( *data, line );
    if ( data->eof( ) ) break;
    Tokenize( line, tokens );
    ForceAssert( tokens.isize( ) > COL );
    ForceAssert( tokens[COL].IsInt( ) );
    lens.push_back( tokens[COL].Int( ) );
  }
  if ( DATA != "-" ) delete data;
  
  sort( lens.begin( ), lens.end( ) );

  // Estimate mean and stdev.
  NormalDistribution nd = SafeMeanStdev( lens );
  float sep = nd.mu_;
  float dev = nd.sigma_;
  float radius = dev * MAX_MULT;
  float bin_begin = sep - radius;
  float bin_end = sep + radius;
  
  // Pick bin size.
  float bin_size = 0.0;
  float pow = 1.0;
  while ( 1 ) {
    for (int ii=0; ii<PROGRESSION.isize( ); ii++) {
      int nbins = int( ( bin_end - bin_begin ) / ( pow * PROGRESSION[ii] ) );
      if ( MIN_BINS <= nbins && nbins < MAX_BINS ) {
	bin_size = pow * PROGRESSION[ii];
	break;
      }
    }
    if ( bin_size > 0 ) break;
    pow *= 10;
  }
    
  // Generate x_bars.
  vec<float> x_bars;
  x_bars.reserve( radius / bin_size );
  int x_bar = 0;
  while ( x_bar < sep - radius )
    x_bar += bin_size;
  while ( x_bar <= sep + radius ) {
    x_bars.push_back( x_bar );
    x_bar += bin_size;
  }
  
  // Build histogram.
  histogram<float> histo;
  for (int ii=0; ii<x_bars.isize( ); ii++)
    histo.AddBin( x_bars[ii] );
  for (int ii=0; ii<lens.isize( ); ii++)
    histo.AddDatum( lens[ii] );
  
  // Histogram's legend.
  color col = black;
  String data_brief;
  if ( DATA != "-" ) {
    data_brief = DATA;
    while ( data_brief.Contains( "/" ) )
      data_brief = data_brief.After( "/" );
  }
  else data_brief = "- (stdin)";
  data_brief = "col " + ToString( COL ) + " from " + data_brief;

  String str_name = "Histogram of " + data_brief;
  String str_mean = ToString( sep ) + " +/- " + ToString( dev );
  String str_sample = "(sample size: " + ToString( lens.isize() ) + " inserts)";
  vec<ns_psplot::freetext> labels( 3 );
  labels[0] = ns_psplot::freetext( str_name, col, 12 );
  labels[1] = ns_psplot::freetext( str_mean, col, 10 );
  labels[2] = ns_psplot::freetext( str_sample, col, 10 );
  
  // Generate eps file.
  ofstream eout( eps_file.c_str( ) );
  histo.PrintAsEps( eout, labels );
  eout.close( );
  
  // Generate text histogram file.
  ofstream tout( histo_file.c_str( ) );
  histo.PrintAsColumns( tout, true, true );
  tout.close( );
  
}
