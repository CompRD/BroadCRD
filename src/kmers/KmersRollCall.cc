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
#include "VecUtilities.h"
#include "kmers/naif_kmer/NaifKmerizer.h" 
#include "kmers/naif_kmer/KernelKmerStorer.h"

typedef Kmer29 Kmer_t;
typedef KmerKmerBiFreq<Kmer_t> KmerRec_t;

/**
 * KmersRollCall
 *
 * Print table of multiplicity of kmers present in reads, reference,
 * or both. Output sent to cout.
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( READS );
  CommandArgument_String( REF );
  CommandArgument_Int_OrDefault( K, 25 );
  CommandArgument_Int_OrDefault( NUM_THREADS, 0 );
  EndCommandArguments;

  NUM_THREADS = configNumThreads( NUM_THREADS );

  // Load.
  vecbvec refreads( REF );
  size_t n_ref = refreads.size( );
  refreads.ReadAll( READS, true );

  // Generate records.
  vec<KmerRec_t> records;
  {
    const bool verbose = false;
    KernelKmerBiStorer<KmerRec_t> storer( refreads, n_ref, K, &records );
    naif_kmerize( &storer, NUM_THREADS, verbose );
  }
  
  // Bins and printable data rows.
  vec<longlong> bins = MkVec( 1l, 2l, 3l, 5l, 10l, 20l, 30l, 50l, 100l );
  vec<longlong> data_row( 3, 0 );   // ref_only, reads_only, both
  vec< vec<longlong> > data( bins.size( ), data_row );
  
  for (size_t ii=0; ii<records.size( ); ii++) {
    const KmerRec_t &record = records[ii];
    vec<unsigned int> acgt = acgt_content( record );

    longlong freq_ref = record.freq_A( );
    longlong freq_reads = record.freq_B( );
    for (int ii=0; ii<bins.isize( ); ii++) {
      longlong threshold = bins[ii];
      if ( freq_ref < threshold && freq_reads < threshold ) break;
      if ( freq_ref >= threshold && freq_reads >= threshold ) data[ii][2]++;
      else if ( freq_reads >= threshold ) data[ii][1]++;
      else data[ii][0]++;
    }
  }

  // Turn into table.
  vec< vec<String> > table;
  vec<String> line;

  line = MkVec( String( "frequency (>=)" ),
		String( "ref_only" ),
		String( "read_only" ),
		String( "both" ),
		String( "total" ) );
  table.push_back( line );

  for (int ii=0; ii<bins.isize( ); ii++) {
    line = MkVec( ToString( bins[ii] ),
		  ToStringAddCommas( data[ii][0] ),
		  ToStringAddCommas( data[ii][1] ),
		  ToStringAddCommas( data[ii][2] ),
		  ToStringAddCommas( BigSum( data[ii] ) ) );
    table.push_back( line );
  }

  // Print.
  cout << "TABLE OF MULTIPLICITY OF " << K << "-MERS\n\n";
  PrintTabular( cout, table, 3, "rrrrr" );
  cout << endl;

}
