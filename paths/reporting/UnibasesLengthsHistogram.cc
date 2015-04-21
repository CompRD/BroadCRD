///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"

#include "Basevector.h"
#include "Histogram.h"

/**
 * UnibasesLengthsHistogram
 *
 * Plot the histogram of unibases (or any fastb file) lengths .
 *
 * OUT_DIR: if empty, it defauls to <UNIBASES>.LengthsHistogram
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( UNIBASES );
  CommandArgument_String_OrDefault( OUT_DIR, "" );
  CommandArgument_Int_OrDefault( LOWEST_BIN, 90 );
  CommandArgument_Int_OrDefault( BIN_SIZE, 1 );
  CommandArgument_Int_OrDefault( N_BINS, 500 );
  EndCommandArguments;

  // Dir and file names.
  String out_dir = OUT_DIR == "" ? UNIBASES + ".LengthsHistogram" : OUT_DIR;

  String log_file = out_dir + "/main.log";
  String eps_file = out_dir + "/histo.eps";
  String txt_file = out_dir + "/histo.txt";
  
  // Output dir, log stream.
  Mkpath( out_dir );
  ofstream log( log_file.c_str( ) );
  PrintCommandPretty( log );
  
  cout << "Sending log to " << log_file << "\n" << endl;

  // Load.
  log << Date( ) << ": loading unibases" << endl;
  vecbvec bases( UNIBASES );
  
  // Bins.
  histogram<int> histo;
  histo.AddLinearProgressionOfBins( LOWEST_BIN, BIN_SIZE, N_BINS );

  // Generate histogram.
  for (size_t ii=0; ii<bases.size( ); ii++) {
    int len = bases[ii].size( );
    histo.AddDatum( len );
  }
  
  // Histogram's legend.
  color col = black;
  String str1 = "Histogram of unibases lengths";
  String str2 = ToString( bases.size( ) ) + " unibases in input";
  
  vec<ns_psplot::freetext> labels( 2 );
  labels[0] = ns_psplot::freetext( str1, col, 12 );
  labels[1] = ns_psplot::freetext( str2, col, 10 );
  
  // Generate eps file.
  ofstream eout( eps_file.c_str( ) );
  histo.PrintAsEps( eout, labels );
  eout.close( );

  // Generate text histogram file.
  ofstream tout( txt_file.c_str( ) );
  histo.PrintAsColumns( tout, true, true );
  tout.close( );
  
  // Done.
  String str_done = Date( ) + ": done";
  cout << str_done << endl;
  log << str_done << endl;

}

