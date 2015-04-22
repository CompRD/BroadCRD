///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Histogram.h"
#include "random/NormalRandom.h"
#include "random/Random.h"

double SimulND( const pair<double,double> &nd );

/**
 * SimulateIndexedMates
 *
 * Simulate relative frequency of amplification from the two end reads
 * of a long molecule. The simulation steps is modeled as follows:
 *
 * 1) molecules are linearly amplified, according to a normal density
 *    function of mean LA_MEAN and deviation LA_DEV. An amplified
 *    molecule is always expected to be shorter than the initial
 *    molecule, and a toss of a coin is performed to check if the
 *    amplification starts from the left or from the right.
 *
 * 2) each amplified molecule therefore contains either the left bit
 *    of an original molecule (and something in the middle of the
 *    original molecule, or the right bit of an original moleculy
 *    (and, again, something in the middle). All of these are PCR
 *    amplified, according to a normal density of mean PCR_MEAN and
 *    deviation PCR_DEV.
 *
 * Output files are saved in OUT_DIR.
 */
int main( int argc, char* argv[] )
{
  RunTime( );

  BeginCommandArguments;

  // Density modeling the linear amplification step.
  CommandArgument_Int( LA_MEAN );
  CommandArgument_Int( LA_DEV );

  // Density modeling the PCR amplification step.
  CommandArgument_Int( PCR_MEAN );
  CommandArgument_Int( PCR_DEV );
  
  // Where output files are saved.
  CommandArgument_String( OUT_DIR );

  // Random seed (used by randomx).
  CommandArgument_Int_OrDefault( SEED, 666 );

  // How many trials to generate, and a dotter to report progress.
  CommandArgument_Int_OrDefault( NSAMPLES, 1000000 );
  CommandArgument_Int_OrDefault( DOTTER, 100000 );
  EndCommandArguments;

  // Dir and file names.
  String log_file = OUT_DIR + "/main.log";
  String eps_file = OUT_DIR + "/histo.eps";
  String txt_file = OUT_DIR + "/histo.txt";
  String raw_file = OUT_DIR + "/trials.output";

  // Random seed.
  srand( SEED );

  // Make outdir, and open log stream.
  Mkpath( OUT_DIR );
  ofstream log( log_file.c_str( ) );
  PrintCommandPretty( log );
  cout << "Sending log to " << log_file << endl;
  
  // Normal distributions (as doubles).
  pair<double,double> nlin = make_pair( double( LA_MEAN ), double( LA_DEV ) );
  pair<double,double> npcr = make_pair( double( PCR_MEAN ), double( PCR_DEV ) );
  
  // The trials.
  triple<double,int,int> empty_trial( 0.0, 0, 0 );
  vec< triple<double,int,int> > trials( NSAMPLES, empty_trial );
  
  // Loop over the number of samples.
  log << Date( ) << ": generating "
      << ToStringAddCommas( NSAMPLES ) << " random trials (.="
      << ToStringAddCommas( DOTTER ) << " trials)" << endl;
  for (int tid=0; tid<NSAMPLES; tid++) {
    if ( tid % DOTTER == 0 ) Dot( log, tid / DOTTER );

    // Linear amplifications.
    int n_lin = SimulND( nlin );
    
    // Loop over all linear amplifications.
    for (int lin_id=0; lin_id<n_lin; lin_id++) {

      // Toss a coin.
      bool left = ( randomx( ) % 2 == 0 );
      int nsequenced = SimulND( npcr );
      
      // Add to proper total.
      if ( left ) ( trials[tid].second ) += nsequenced;
      else ( trials[tid].third ) += nsequenced;
    }
    
    // Make sure left <= right.
    if ( trials[tid].third < trials[tid].second )
      swap( trials[tid].second, trials[tid].third );
    
    // Compute ratio.
    trials[tid].first 
      = ( trials[tid].third < 1 )
      ? 0.0
      : (double)trials[tid].second / (double)trials[tid].third;
    
  }
  log << endl;
  
  // Sort.
  log << Date( ) << ": sorting trials" << endl;
  sort( trials.begin( ), trials.end( ) );

  // Dump trials to file.
  {
    ofstream out( raw_file.c_str( ) );
    for (int ii=0; ii<trials.isize( ); ii++)
      out << trials[ii].first << "\t"
       	  << trials[ii].second << "\t"
       	  << trials[ii].third << "\n";
    out.close( );
  }
  
  // Build histogram (note last bin is 0.99, not 1.0).
  vec<float> bins;
  for (int ii=0; ii<100; ii++)
    bins.push_back( float( ii ) / 100.0 );

  histogram<float> histo;
  histo.AddBins( bins );
  for (int ii=0; ii<trials.isize( ); ii++)
    histo.AddDatum( trials[ii].first );
    
  // Plot histogram.
  using namespace ns_psplot;

  const String str_nsamples = ToStringAddCommas( NSAMPLES );
  const String str_lin = ToString( LA_MEAN ) + " +/- " + ToString( LA_DEV );
  const String str_pcr = ToString( PCR_MEAN ) + " +/- " + ToString( PCR_DEV );

  vec<ns_psplot::freetext> labels( 4 );
  labels[0] = freetext( "Histogram of simulated trials", black, 16 );
  labels[1] = freetext( str_nsamples + " random samples", black, 14 );
  labels[2] = freetext( "linear amplification:  " + str_lin, black, 14 );
  labels[3] = freetext( "pcr amplification:  " + str_pcr, black, 14 );
  
  ofstream eps_out( eps_file.c_str( ) );
  histo.PrintAsEps( eps_out, labels, 2 );
  eps_out.close( );  

  ofstream txt_out( txt_file.c_str( ) );
  txt_out << "Histogram of simulated trials\n"
	  << str_nsamples << " random samples\n"
	  << "linear amplification:  " << str_lin << "\n"
	  << "pcr amplification:  " << str_pcr << "\n"
	  << endl;

  histo.PrintAsColumns( txt_out );
  txt_out.close( );

  // Done.
  String date = Date( );
  cout << date << ": SimulateIndexedMates done" << endl;
  log << date << ": SimulateIndexedMates done" << endl;
  log.close( );

}

/**
 * SimulND
 *
 * Just a wrapper around FastNormal
 */
double SimulND( const pair<double,double> &nd )
{
  return nd.first + nd.second * FastNormal( );
}

