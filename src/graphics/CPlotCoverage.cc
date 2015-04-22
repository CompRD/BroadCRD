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
#include "graphics/BarGraph.h"
#include "graphics/CContigPlotter.h"
#include "graphics/CPlotCoverage.h"
#include "reporting/SnapToGrid.h"
#include "util/PickTic.h"

/**
 * CPlotCoverage
 * Constructor
 */
CPlotCoverage::CPlotCoverage( int win,
			      float cap,
			      const String &name,
			      const String &logfile ) :
  win_size_ ( win ),
  hop_size_ ( win / 2 ),
  cap_ ( cap ),
  mult_ ( 3.0 ),
  name_ ( name ),
  mode_ ( -1 )
{
  ForceAssertGe( hop_size_, 1 );
  ForceAssertGe( win_size_, 1 );

  log_.open( logfile.c_str( ) );
}

/**
 * CPlotCoverage
 * Destructor
 */
CPlotCoverage::~CPlotCoverage( )
{
  log_.close( );
}

/**
 * CPlotCoverage
 * SetFromSeqIntervals
 *
 * The input consists of a vector of seq_intervals generated, for instance,
 * by CoverageAnalyzer::GetAllCoverages( ).
 */
void CPlotCoverage::SetFromSeqIntervals( const vec<int> &tlens,
					 const vec<int> &fcov,
					 const vec<seq_interval> &coverages )
{
  // Clean and resize.
  this->ResizeAverages( tlens );

  // Loop over all targets.
  for (int target_id=0; target_id<(int)tlens.size( ); target_id++) {
    int target_len = tlens[target_id];

    // Test for early exit.
    int firstcov = fcov[ target_id ];
    if ( firstcov < 0 || target_len < win_size_ )
      continue;

    // Fill cov_.
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

    // File averages_.
    this->GenerateAverages( target_id, raw );
  }
}

/**
 * CPlotCoverage
 * SetFromVecvec
 *
 * The input consists of a vecqualvector representig (capped) coverage,
 * as from the output of GenerateCoverageVecvec.
 */
void CPlotCoverage::SetFromVecvec( const vecqualvector &vcov )
{
  // Clean and resize.
  vec<int> tlens( vcov.size( ) );
  for (vecqvec::size_type ii=0; ii<vcov.size( ); ii++)
    tlens[ii] = vcov[ii].size( );
  this->ResizeAverages( tlens );

  // Loop over all targets.
  for (int target_id=0; target_id<(int)tlens.size( ); target_id++) {
    int target_len = vcov[target_id].size( );
    if ( target_len < win_size_ )
      continue;

    // Fill cov_.
    vec<float> raw( target_len, 0 );
    for (int ii=0; ii<target_len; ii++)
      raw[ii] = (float)(int)vcov[target_id][ii];

    // File averages_.
    this->GenerateAverages( target_id, raw );
  }
}

/**
 * CPlotCoverage
 * PlotTargets
 */
void CPlotCoverage::PlotTargets( const String &out_dir, bool histo_only )
{
  ForceAssert( IsDirectory( out_dir ) );

  longlong tot_windows = 0;
  longlong tot_tagged = 0;

  String histo_file = out_dir + "/histogram.eps";
  float mode = this->GlobalMode( histo_file );
  float max_aver = mult_ * mode;

  if ( histo_only )
    return;

  // Loop over all targets.
  for (int target_id=0; target_id<(int)averages_.size( ); target_id++) {
    int length = tlens_[target_id];
    const vec<float> &aver = averages_[target_id];
    const String str_target = ToString( target_id );
    const String outfile = out_dir + "/target_" + str_target + ".eps";

    longlong tagged = 0;

    log_ << "Plotting target " << target_id
	 << " (len=" << length
	 << ", n_points=" << aver.size( )
	 << ")\n";

    // Data plot.
    CDataPlotter data( red );
    data.Reserve( aver.size( ) );

    float data_max = 0;
    for (int ii=0; ii<(int)aver.size( ); ii++) {
      float xpt = ii * hop_size_;
      float ypt = aver[ii];

      if ( ypt > max_aver ) {
	tagged++;
	int begin = ii * hop_size_;
	int end = Min( begin + win_size_, length );
	float excess = aver[ii] / mode;

	log_ << " t_" << target_id
	     << " [" << begin / 1000
	     << " Kb, " << end / 1000
	     << " Kb) cov = " << ToString( aver[ii], 2 )
	     << "X ie " << ToString( excess, 2 )
	     << " times the mode\n";
      }

      data.AddPoint( xpt, ypt );
      data_max = Max( data_max, ypt );
    }

    tot_windows += aver.size( );
    tot_tagged += tagged;
    log_ << "Stats for target_" << target_id
	 << ": " << aver.size( )
	 << " data points, of which " << tagged
	 << " where tagged as high coverage\n" << endl;

    // Mode plots.
    CDataPlotter mplot( blue );
    mplot.AddPoint( 0, mode );
    mplot.AddPoint( (float)length, mode );

    CDataPlotter mhplot( blue );
    mhplot.AddPoint( 0, max_aver );
    mhplot.AddPoint( (float)length, max_aver );

    // Tics
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

    // Plot.
    CContigPlotter plotter( length );
    plotter.AddData( data );
    plotter.AddData( mplot );
    plotter.AddData( mhplot );
    plotter.AddHorizontalTics( xtics );
    plotter.AddVerticalTics( ytics );

    plotter.GoFigure( outfile, true );

  } // Loop over targets.

  // Final stats.
  double ratio = 100.0 * SafeQuotient( tot_tagged, tot_windows );

  log_ << "Final statistics:\n"
       << " window size: " << win_size_ << "\n"
       << " slide size (hop amount): " << hop_size_ << "\n"
       << " mode multiplier to tag window: " << ToString( mult_ ) << "\n"
       << " number of data points total: " << tot_windows << "\n"
       << " number of windows tagged as high coverage: " << tot_tagged << "\n"
       << " percent of windows tagged: " << ToString( ratio, 2 ) << "%\n"
       << "\n" << endl;

}

/**
 * CPlotCoverage
 * GolbalMode
 *
 * If coverage is very low, many windows will end up with zero (or almost
 * zero) coverage. The distribution is therefore generated only on the
 * windows for which coverage is bigger than a fixed given epsilon.
 */
float CPlotCoverage::GlobalMode( const String &histo_file )
{
  double epsilon = 0.01;

  if ( mode_ > -1 )
    return mode_;

  longlong nobj = 0;
  for (int ii=0; ii<(int)averages_.size( ); ii++)
    nobj += averages_[ii].size( );

  vec<float> aver;
  aver.reserve( nobj );
  for (int ii=0; ii<(int)averages_.size( ); ii++)
    for (int jj=0; jj<(int)averages_[ii].size( ); jj++)
      if ( averages_[ii][jj] >= epsilon )
	aver.push_back( averages_[ii][jj] );

  NormalDistribution nd = SafeMeanStdev( aver );
  float mean = nd.mu_;
  float stdev = nd.sigma_;

  // Generate fine-grid histogram to estimate mode.
  float local_cap = mean + 12.0 * stdev;

  vec<float> xbins;
  float bin = 0;
  float increm = 0.01;
  while ( bin * increm < local_cap ) {
    xbins.push_back( bin * increm );
    bin++;
  }

  vec<float> ybins( xbins.size( ), 0 );
  for (int ii=0; ii<(int)aver.size( ); ii++) {
    int pos = PlaceInBin( aver[ii], xbins );
    ybins[pos] += 1;
  }

  mode_ = 0;
  float maxbin = ybins[0];
  for (int ii=0; ii<(int)ybins.size( ); ii++) {
    if ( ybins[ii] > maxbin ) {
      maxbin = ybins[ii];
      mode_ = xbins[ii];
    }
  }
  String str_mode = ToString( mode_, 2 );

  // Generate histogram for plotting.
  xbins.clear( );
  bin = 0;
  increm = 0.1;
  while ( bin * increm <= cap_ + epsilon ) {
    xbins.push_back( bin * increm );
    bin++;
  }

  ybins.clear( );
  ybins.resize( xbins.size( ), 0 );
  for (int ii=0; ii<(int)aver.size( ); ii++) {
    int pos = PlaceInBin( aver[ii], xbins );
    ybins[pos] += 1;
  }

  log_ << "\nHistogram of coverages (see also ./" << histo_file << "):\n";
  for (int ii=0; ii<(int)xbins.size( ); ii++)
    log_ << " " << xbins[ii] << "\t" << ybins[ii] << "\n";
  log_ << "\n";

  String str_eps = ToString( epsilon, 2 );
  String str_mean = ToString( mean ) + " +/- " + ToString( stdev );
  log_ << "Distrib. is generated on windows with cov >= " << str_eps << "\n"
       << " data points in input: " << nobj << "\n"
       << " data points used: " << aver.size( ) << "\n"
       << " estimated mean: " << str_mean << "\n"
       << " estimated mode: " << str_mode << "\n"
       << endl;

  float low = 0;
  float high = Max( ybins );

  float max_bar_height = 150;
  float bar_width = 6;
  float bar_sep = 1.5;
  short font_size = 7;
  color col = gray;

  String decor = "Histogram of coverages for " + name_;
  String str_win = ToString( win_size_ );

  vec<ns_psplot::freetext> d_labels( 4 );
  d_labels[0] = ns_psplot::freetext( decor, black, 12 );
  d_labels[1] = ns_psplot::freetext( "window size: " + str_win, black, 9 );
  d_labels[2] = ns_psplot::freetext( "mean: " + str_mean, black, 9 );
  d_labels[3] = ns_psplot::freetext( "mode: " + str_mode, black, 9 );

  ns_psplot::bargraph d_graph = BarPlot( ybins, low, high, bar_width, bar_sep,
					 col, font_size, max_bar_height,
					 d_labels, xbins, true);

  ofstream out( histo_file.c_str( ) );
  out << d_graph << "\n";
  out.close( );

  // Ok, return.
  return mode_;
}

/**
 * CPlotCoverage
 * ResizeAverages
 */
void CPlotCoverage::ResizeAverages( const vec<int> &tlens )
{
  mode_ = -1;
  tlens_ = tlens;
  averages_.clear( );
  averages_.resize( tlens.size( ) );

  for (int target_id=0; target_id<(int)tlens.size( ); target_id++) {
    int data_size = 0;
    while ( data_size * hop_size_ < (int)averages_[target_id].size( ) )
      data_size++;
    averages_[target_id].reserve( data_size );
  }
}

/**
 * CPlotCoverage
 * GenerateAverages
 */
void CPlotCoverage::GenerateAverages( int target_id, const vec<float> &raw )
{
  int post = 0;
  while ( post * hop_size_ < tlens_[target_id] ) {
    longlong totcov = 0;
    int begin = post * hop_size_;
    int end = Min( begin + win_size_, tlens_[target_id] );
    for (int ii=begin; ii<end; ii++)
      totcov += (int)raw[ii];
    float average = (float)totcov / (float)( end - begin );
    averages_[target_id].push_back( average );
    post++;
  }

  NormalDistribution nd = SafeMeanStdev( averages_[target_id] );
  String str_mean = ToString( nd.mu_ ) + " +/- " + ToString( nd.sigma_ );

  log_ << " mean coverage for target_" << target_id
       << ": " << str_mean
       << " (over " << averages_[target_id].size( )
       << " windows)\n" << flush;
}

