///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Basevector.h"
#include "Bitvector.h"
#include "FeudalMimic.h"
#include "Qualvector.h"
#include "VecUtilities.h"
#include "feudal/BinaryStream.h"
#include "kmers/KmersTaggerCore.h"
#include "kmers/ReadPather.h"
#include "kmers/ReadPatherDefs.h"
#include "kmers/naif_kmer/NaifKmerizer.h" 
#include "kmers/naif_kmer/KernelKmerStorer.h"
#include "util/RunCommand.h"

void PrintTable( const vec< vec<size_t> > &raw,
		 const vec<int> &min_freqs,
		 const String &title,
		 ostream &out );

/**
 * EvaluateKmersTagger
 *
 * Run KmerstTaggerCore on a set of reads, and generate coverage of
 * references based on various MIN_FREQ sizes. Output sent to cout
 * (but several files are cached in OUT_DIR).
 *
 * REMARKS:
 * 
 * 1. The kmer size for error correction is hard coded (the argument K
 *    is used by the evaluation code only).
 *
 * 2. KmersTaggerCore uses Ted's ReadPather to tag kmers and remove
 *    untrusted portions of reads, while this evaluation code uses
 *    Filipe's NaifKmerizer.
 *
 * REF: head of reference
 * READS: head of reads
 * OUT_DIR: where cached files are saved
 * KEVAL: kmer size (for evaluation, the error correction K is hard coded)
 * MIN_READLEN: do not consider short segments (kmer length)
 * MIN_FREQS: parsed with ParseIntSet
 * FORCE: do not used cached data.
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments; 
  CommandArgument_String( REF );
  CommandArgument_String( READS );
  CommandArgument_String_OrDefault( OUT_DIR, READS + ".EvaluateKmersTagger" );
  CommandArgument_Int_OrDefault( KEVAL, 25 );
  CommandArgument_Int_OrDefault( MIN_READLEN, 12 );
  CommandArgument_String_OrDefault( MIN_FREQS, "{1,2,3,4,5}");
  CommandArgument_UnsignedInt_OrDefault( NUM_THREADS, 0 );
  CommandArgument_Bool_OrDefault( FORCE, True );
  CommandArgument_Bool_OrDefault( VALIDATE, False );
  EndCommandArguments;  

  // Constants.
  const unsigned K = 21;

  // Thread control.
  NUM_THREADS = configNumThreads( NUM_THREADS );

  // File names.
  String strK = ToString( K );
  String ref_file = REF + ".fastb";
  String refa_file = REF + ".fasta";
  String refnames_file = REF + ".names";
  String reads_file = READS + ".fastb";
  String dict_file = READS + ".k" + strK + ".ug.dict";

  String log_file = OUT_DIR + "/main.log";
  
  vec<String> needed;
  needed.push_back( ref_file );
  needed.push_back( refa_file );
  needed.push_back( reads_file );

  Mkpath( OUT_DIR );
  
  // Load.
  vec<int> min_freqs;
  ParseIntSet( MIN_FREQS, min_freqs );

  vec<String> ref_names;
  if ( ! IsRegularFile( refnames_file ) ) {
    int n_refs = MastervecFileObjectCount( ref_file );
    ref_names.reserve( n_refs );
    String line;
    ifstream in( refa_file.c_str( ) );
    while( in ) {
      getline( in, line );
      if ( ! in ) break;
      if ( ! line.Contains( ">" ) ) continue;
      String name = line.After( ">" );
      while ( name.Contains( " " ) )
	name = name.Before( " " );
      ref_names.push_back( name );
    }
    in.close( );
    WRITE( refnames_file, ref_names );
  }
  else
    READX( refnames_file, ref_names );
  
  // Run KmersTaggerCore on all min_freqs.
  bool run_kmers_tagger = false;
  for (size_t freq_id=0; freq_id<min_freqs.size( ); freq_id++) {
    const String str_minF = ToString( min_freqs[freq_id] );
    const String out_head = OUT_DIR + "/segs_MF" + str_minF;
    const String segs_MF_file = out_head + ".fastb";
    if ( FORCE || ! IsRegularFile( segs_MF_file ) ) {
      run_kmers_tagger = true;
      break;
    }
  }
  
  if ( run_kmers_tagger ) { 
    cout << Date( ) << ": loading reads" << endl;
    VirtualMasterVec<bvec> reads( reads_file );
    
    cout << Date( ) << ": loading reference" << endl;
    VirtualMasterVec<bvec> ref( ref_file );
    
    KmerDict<K> dictW( 0 );
    KmerDict<K> dictR( 0 );
    bool write_dict = ( FORCE || ! IsRegularFile( dict_file ) );
    if ( write_dict ) {
      cout << Date( ) << ": generating (and saving) dict" << endl;
      dictW.process( reads, VALIDATE, NUM_THREADS );
      BinaryWriter::writeFile( dict_file, dictW );
    }
    else {
      cout << Date( ) << ": reading dict" << endl;
      BinaryReader::readFile( dict_file, &dictR );
    }
    const KmerDict<K> &dict = write_dict ? dictW : dictR;
    cout << endl;
    
    // Loop over all min_freqs.
    for (size_t freq_id=0; freq_id<min_freqs.size( ); freq_id++) {
      const int minF = min_freqs[freq_id];
      const int minRL = MIN_READLEN;
      const bool verb = false;
      
      // Local file names.
      const String str_minF = ToString( minF );
      const String out_head = OUT_DIR + "/segs_MF" + str_minF;
      const String stats_MF_file = out_head + ".stats";
      const String segs_MF_file = out_head + ".fastb";
      
      cout << Date( ) << ": starting MIN_FREQ = " << str_minF << endl;
      
      // Generate valid segments.
      if ( FORCE || ! IsRegularFile( segs_MF_file ) ) {
	if ( minF == 1 ) {
	  cout << Date( ) << ": copying reads" << endl;
	  Cp( reads_file, segs_MF_file );
	}
	else {
	  cout << Date( ) << ": running KmersTaggerCore" << endl;
	  vecbvec segs;
	  KmersTaggerCore<K> ( minF, minRL, verb, dict, reads, segs, cout );
	  cout << Date( ) << ": saving valid segments" << endl;
	  segs.WriteAll( segs_MF_file );
	}
      }
      cout << endl;
    }
  }
  
  // Load reference.
  cout << Date( ) << ": loading reference\n" << endl;
  vecbvec ref( ref_file );
  size_t n_ref = ref.size( );
    
  // Generate stats files on all min_freqs.
  bool run_naif = false;
  for (size_t freq_id=0; freq_id<min_freqs.size( ); freq_id++) {
    const String str_minF = ToString( min_freqs[freq_id] );
    const String out_head = OUT_DIR + "/segs_MF" + str_minF;
    const String stats_MF_file = out_head + ".stats";
    if ( FORCE || ! IsRegularFile( stats_MF_file ) ) {
      run_naif = true;
      break;
    }
  }
  
  if ( run_naif ) {
    for (size_t freq_id=0; freq_id<min_freqs.size( ); freq_id++) {
      const int minF = min_freqs[freq_id];
      const bool verb = false;
      
      // Local file names.
      const String str_minF = ToString( minF );
      const String out_head = OUT_DIR + "/segs_MF" + str_minF;
      const String segs_MF_file = out_head + ".fastb";
      const String stats_MF_file = out_head + ".stats";
      
      // Run the Naif Kmerizer package.
      if ( FORCE || ! IsRegularFile( stats_MF_file ) ) {
	typedef Kmer29 Kmer_t;
	typedef KmerKmerBiFreq<Kmer_t> KmerRec_t;
	const bool verbose = false;

	cout << Date( ) << ": loading segs_MF" + str_minF << endl;
	vecbvec alldata = ref;
	alldata.ReadAll( segs_MF_file, true );

	cout << Date( ) << ": naif kmerizing MIN_FREQ = " << minF << endl;
	vec<KmerRec_t> records;
	KernelKmerBiStorer<KmerRec_t> storer( alldata, n_ref, KEVAL, &records );
	naif_kmerize( &storer, NUM_THREADS, verbose );
	
	size_t in_both = 0;
	size_t ref_only = 0;
	size_t reads_only = 0;
	for (size_t ii=0; ii<records.size( ); ii++) {
	  const KmerRec_t &record = records[ii];
	  longlong freq_ref = record.freq_A( );
	  longlong freq_reads = record.freq_B( );
	  ForceAssert( freq_ref > 0 || freq_reads > 0 );
	  if ( freq_reads > 0 && freq_ref > 0 ) in_both++;
	  else if ( freq_reads > 0 ) reads_only++;
	  else ref_only++;
	}
	
	vec< vec<String> > table;
	vec<String> line = MkVec( String( "#" ),
				  String( "in_both" ),
				  String( "reads_only" ),
				  String( "ref_only" ) );
	table.push_back( line );
	
	line = MkVec( String( "" ),
		      ToString( in_both ),
		      ToString( reads_only ),
		      ToString( ref_only ) );
	table.push_back( line );
	
	String rjust = "rrrr";
	ofstream out( stats_MF_file.c_str( ) );
	PrintTabular( out, table, 3, rjust );
	out.close( );

	cout << endl;
      }
    }
  }

  // raw[freq_id][0]: in both, [1]: reads_only, [2]: ref_only
  vec< vec<size_t> > raw;
  raw.resize( min_freqs.size( ) );
  for (size_t ii=0; ii<min_freqs.size( ); ii++)
    raw[ii].resize( 3, 0 );
  
  // Parse stats files.
  for (size_t ii=0; ii<min_freqs.size( ); ii++) {
    String line;
    String str_mf = ToString( min_freqs[ii] );
    String stats_MF_file = OUT_DIR + "/segs_MF" + str_mf + ".stats";
    
    ifstream in( stats_MF_file.c_str( ) );
    while ( in ) {
      getline( in, line );
      in >> raw[ii][0] >> raw[ii][1] >> raw[ii][2];
    }
    in.close( );
  }
  
  // Generate table. 
  String title = "TABLE OF GENOMIC VS NON-GENOMIC KMERS";
  PrintTable( raw, min_freqs, title, cout );

  // Done.
  cout << Date( ) << ": EvaluateKmersTagger done" << endl;

}



/**
 * PrintTable
 */
void PrintTable( const vec< vec<size_t> > &raw,
		 const vec<int> &min_freqs,
		 const String &title,
		 ostream &out )
{
  vec< vec<String> > table;
  
  String rjust = "rr";
  vec<String> line;
  line.push_back( String( "" ) );
  line.push_back( String( "reference" ) );
  for (size_t ii=0; ii<min_freqs.size( ); ii++) {
    line.push_back( "mf_" + ToString( min_freqs[ii] ) );
    rjust += 'r';
  }
  table.push_back( line );
  
  {
    line.clear( );
    line.push_back( "%% genomic missing" );
    line.push_back( "0.00" );
    for (size_t ii=0; ii<min_freqs.size( ); ii++) {
      size_t tot = raw[ii][0] + raw[ii][2];
      size_t found = raw[ii][0];
      double ratio = 10000. * ( 1. - SafeQuotient( found, tot ) );
      line.push_back( ToString( ratio, 2 ) );
    }
    table.push_back( line );
  }

  {
    line.clear( );
    line.push_back( "genomic missing" );
    line.push_back( "0" );
    for (size_t ii=0; ii<min_freqs.size( ); ii++)
      line.push_back( ToStringAddCommas( raw[ii][2] ) );
    table.push_back( line );
  }
  
  {
    line.clear( );
    line.push_back( "genomic found" );
    line.push_back( ToStringAddCommas( raw[0][0] + raw[0][2] ) );
    for (size_t ii=0; ii<min_freqs.size( ); ii++)
      line.push_back( ToStringAddCommas( raw[ii][0] ) );
    table.push_back( line );
  }
  
  {
    line.clear( );
    line.push_back( "non genomic" );
    line.push_back( "0" );
    for (size_t ii=0; ii<min_freqs.size( ); ii++)
      line.push_back( ToStringAddCommas( raw[ii][1] ) );
    table.push_back( line );
  }

  out << title << "\n\n";
  PrintTabular( out, table, 3, rjust );
  out << endl;
  
}
