/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "Histogram.h"
#include "VecUtilities.h"
#include "feudal/BinaryStream.h"
#include "kmers/ReadPather.h"
#include "kmers/ReadPatherDefs.h"

/**
 * KmersHistogram
 *
 * Use Ted's pathing code ("new", as to May 2012), to generate an
 * histogram of kmers' frequency.  This code is here for two reasons:
 * first, this is the most efficient way I know of to generate such an
 * histogram; and secondly, as an example of how to use Ted's code.
 *
 * HEAD_IN: it needs the .fastb and .k<K>.ug.dict (<K> hard coded below)
 * HEAD_OUT: it generates both txt and eps output
 * NUM_THREADS: use all if 0
 * FORCE: do not used cached dictionary
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( HEAD_IN );
  CommandArgument_String_OrDefault( HEAD_OUT, HEAD_IN );
  CommandArgument_Int_OrDefault( NUM_THREADS, 0 );
  CommandArgument_Bool_OrDefault( FORCE, True );
  CommandArgument_Bool_OrDefault( VALIDATE, False );
  EndCommandArguments;
 
  // Thread control.
  NUM_THREADS = configNumThreads( NUM_THREADS );

  // Constants.
  const unsigned K = 24;
  const size_t MAX_BINS = 2000;
  const double PERCENT_TO_SHOW = 99.5;
    
  // Dir and file names.
  String strK = ToString( K );
  String bases_file = HEAD_IN + ".fastb";
  String dict_file = HEAD_IN + ".k" + strK + ".ug.dict";

  String out_txt_file = HEAD_OUT + ".KmersHistogram.txt";
  String out_eps_file = HEAD_OUT + ".KmersHistogram.eps";

  // Load.
  cout << Date( ) << ": loading bases" << endl;
  VirtualMasterVec<bvec> vmv_bases( bases_file );
  
  KmerDict<K> dictW( 0 );
  KmerDict<K> dictR( 0 );
  bool write_dict = ( FORCE || ! IsRegularFile( dict_file ) );
  if ( write_dict ) {
    cout << Date( ) << ": generating (and saving) dict" << endl;
    dictW.process( vmv_bases, VALIDATE, NUM_THREADS );
    BinaryWriter::writeFile( dict_file, dictW );
  }
  else {
    cout << Date( ) << ": reading dict" << endl;
    BinaryReader::readFile( dict_file, &dictR );
  }
  const KmerDict<K> &dict = write_dict ? dictW : dictR;
  const size_t n_kmers = dict.size( );
  const String str_n_kmers = ToStringAddCommas( n_kmers );

  // Parsing reads (copied from assessDictionary in paths/long/Friends.cc).
  cout << Date( ) << ": parsing " << str_n_kmers << " kmers" << endl;
  vec<size_t> k_counts( MAX_BINS + 1 );
  typedef KmerDict<K>::OCItr OCItr;
  typedef KmerDict<K>::ICItr ICItr;
  for (OCItr oItr(dict.begin( )), oEnd(dict.end( )); oItr!=oEnd; ++oItr) {
    for (ICItr itr(oItr->begin( )), end(oItr->end( )); itr!=end; ++itr) {
      KDef const *pKDef = dict.lookup( *itr );
      ForceAssert( pKDef );
      size_t count = Min( pKDef->getCount( ), MAX_BINS );
      k_counts[count] += 1;
    }
  }
  ForceAssertEq( BigSum( k_counts ), n_kmers );

  // Find a reasonable cut for the bins to show.
  cout << Date( ) << ": generating histogram" << endl;
  size_t last_bin = 0;
  size_t seen_kmers = 0;
  for (size_t ii=0; ii<=MAX_BINS; ii++) {
    seen_kmers += k_counts[ii];
    if ( 100. * SafeQuotient( seen_kmers, n_kmers ) >= PERCENT_TO_SHOW ) {
      last_bin = ii + 1 ;
      break;
    }
  }
  while ( last_bin % 10 != 0 ) last_bin++;
  
  // Generate histogram.
  histogram<int> histo;
  histo.AddLinearProgressionOfBins( 0, 1, last_bin + 1 );
  for (size_t ii=0; ii<=MAX_BINS; ii++)
    histo.AddIdenticalData( (int)ii, k_counts[ii] );

  // Save output files.
  ofstream txt_out( out_txt_file.c_str( ) );
  histo.PrintAsColumns( txt_out, true, true, true );
  txt_out.close( );

  typedef ns_psplot::freetext freetxt;
  vec<freetxt> labels;
  labels.push_back( freetxt( "Kmers distribution" , black, 12 ) );
  labels.push_back( freetxt( bases_file, black, 10 ) );
  labels.push_back( freetxt( str_n_kmers + " kmers in input", black, 10 ) );
  ofstream eps_out( out_eps_file.c_str( ) );
  histo.PrintAsEps( eps_out, labels, 0 );
  eps_out.close( );
  
  // Done.
  cout << Date( ) << ": done" << endl;

}

