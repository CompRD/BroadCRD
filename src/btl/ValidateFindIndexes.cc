///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"

#include "Fastavector.h"
#include "Histogram.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "util/RunCommand.h"
#include "paths/reporting/ReftigUtils.h"
#include <omp.h>
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

/**
 * ValidateFindIndexes
 *
 * Load the output from FindIndexes, and validate the genomic chunks
 * against a given reference. Reads will be saved in pairs (where
 * entries 2*i and 1 + 2*i are mates to each other).
 *
 * IN_DIR: where input files are
 * REF_L: lookup table of reference
 */
int main( int argc, char* argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( IN_DIR );
  CommandArgument_String( REF_L );

  // If reference is circular.
  CommandArgument_Bool_OrDefault( CIRCULAR, False );

  // Define range of valid separations, and histogram related args.
  CommandArgument_Int_OrDefault( MIN_SEP, 0 );
  CommandArgument_Int_OrDefault( MAX_SEP, 75000 );
  CommandArgument_Int_OrDefault( BIN_SIZE, 500 );

  // Aligner's args.
  CommandArgument_Int_OrDefault( K, 26 );
  CommandArgument_Int_OrDefault( NUM_THREADS, 0 );
  CommandArgument_Bool_OrDefault( FORCE, True );

  EndCommandArguments;
  
  // Check args.
  if ( K != 26 ) {
    cout << "Fatal error: for now only allowed value for K is 26.\n" << endl;
    return 1;
  }

  // Dir and file names.
  const String reads_head = "indexed_reads";

  String full_reads_head = IN_DIR + "/" + reads_head;
  String rbases_file = full_reads_head + ".fasta";
  String barcounts_file = full_reads_head + ".barcounts";

  String out_dir = full_reads_head + ".Validate";
  String log_file = out_dir + "/main.log";
  String reads_file = out_dir + "/reads.fastb";
  String names_file = out_dir + "/reads.names";
  String hits_file = out_dir + "/reads.qlt";
  String best_hits_file = out_dir + "/select.qlt";
  String seps_file = out_dir + "/select.seps";
  String histo_eps_file = out_dir + "/histogram.eps";
  String histo_txt_file = out_dir + "/histogram.txt";
  
  vec<String> needed;
  needed.push_back( REF_L );
  needed.push_back( rbases_file );
  needed.push_back( barcounts_file );

  if ( ! CheckFilesExist( needed, &cout ) ) return 1;
  
  Mkpath( out_dir );
  
  // Thread control (needed by FastAlignShortReads).
  NUM_THREADS = configNumThreads( NUM_THREADS );
  omp_set_num_threads( NUM_THREADS );

  // Log stream.
  ofstream log( log_file.c_str( ) );
  PrintCommandPretty( log );
  cout << " Sending log to " << log_file << "\n" << endl;
  
  // Load barcounts file.
  log << Date( ) << ": loading barcounts" << endl;
  vec<long> counts;
  vec<String> indexes;
  {
    String name = "";
    int cnt = 0;
    ifstream in( barcounts_file.c_str( ) );
    while( in ) {
      in >> name >> cnt;
      if ( !in ) break;
      
      if ( cnt == 2 ) indexes.push_back( name );
      if ( counts.isize( ) <= cnt ) counts.resize( cnt + 1 );
      counts[cnt] += 1;
    }
    in.close( );
  }

  // Print table of multiplicities.
  vec< vec<String> > table;
  
  long tot = BigSum( counts );
  table.push_back( MkVec( String( "tot" ),
			  ToStringAddCommas( tot ),
			  String( "100.0%" ) ) );
  
  for (int ii=1; ii<counts.isize( ); ii++) {
    double ratio = counts[ii] == 0 ? 0. : SafeQuotient( counts[ii], tot );
    table.push_back( MkVec( ToStringAddCommas( ii ),
			    ToStringAddCommas( counts[ii] ),
			    ToString( 100. * ratio, 1 ) + "%" ) );
  }
  
  log << "\n"
      << "TABLE OF MULTIPLICITIES\n"
      << "\n";
  PrintTabular( log, table, 3, "rrr" );
  log << endl;
  
  // Convert reads to fastb.
  if ( FORCE || ! IsRegularFile( reads_file ) ) {
    log << Date( ) << ": loading genomic chunks" << endl;

    vec<fastavector> fastas;
    vec<String> names;
    LoadFromFastaFile( rbases_file, fastas, names );
    
    log << Date( ) << ": selecting barcodes with exactly two chunks" << endl;
    vecbvec bases;
    vec<String> ids;
    bases.reserve( 2 * indexes.size( ) );
    ids.reserve( 2 * indexes.size( ) );
    sort( indexes.begin( ), indexes.end( ) );
    for (size_t ii=0; ii<fastas.size( ); ii++) {
      String name = names[ii].After( "_" ).Before( "_" );
      if ( binary_search( indexes.begin( ), indexes.end( ), name ) ) {
	bases.push_back( fastas[ii].ToBasevector( ) );
	ids.push_back( names[ii] );
      }
    }
    
    log << Date( ) << ": saving selected chunks" << endl;
    bases.WriteAll( reads_file );
    WRITE( names_file, ids );
  }
  
  // Align reads, and filter aligns.
  if ( FORCE || ! IsRegularFile( best_hits_file ) ) {
    const size_t nreads = MastervecFileObjectCount( reads_file );
    const size_t npairs = nreads / 2;

    ForceAssert( nreads % 2 == 0 );

    log << Date( ) << ": aligning reads" << endl;
    vec<look_align> hits;
    GetAlignsFast( K, reads_file, REF_L, hits_file, hits, !FORCE, out_dir );
    
    log << Date( ) << ": filtering alignments" << endl;
    vec<int> to_hit( nreads, -1);
    for (size_t ii=0; ii<hits.size( ); ii++) {
      int id = hits[ii].query_id;
      if ( to_hit[id] == -1 ) to_hit[id] = ii;
      else to_hit[id] = -2;
    }

    // Keep track of discards
    //   0:   valid
    //   1:   one or both end read unaligned
    //   2:   one or both end reads multiply placed
    //   3:   on different targets
    //   4:   wrong orientation
    //   5:   duplicate (should be rare)
    vec<int> invalids( 6, 0 );
  
    ofstream best_out( best_hits_file.c_str( ) );
    ofstream seps_out( seps_file.c_str( ) );
    for (size_t ii=0; ii<npairs; ii++) {
      int id1 = ii * 2;
      int id2 = id1 + 1;
      
      // One or both end reads unaligned.
      if ( to_hit[id1] == -1 || to_hit[id2] == -1 ) {
	invalids[1] += 1;
	continue;
      }

      // One or both end reads multiply aligned.
      if ( to_hit[id1] == -2 || to_hit[id2] == -2 ) {
	invalids[2] += 1;
	continue;
      }

      // On different targets.
      const look_align &hit1 = hits[ to_hit[id1] ];
      const look_align &hit2 = hits[ to_hit[id2] ];
      if ( hit1.TargetId( ) != hit2.TargetId( ) ) {
	invalids[3] += 1;
	continue;
      }

      // Same orientation.
      if ( hit1.Rc1( ) == hit2.Rc1( ) ) {
	invalids[4] += 1;
	continue;
      }
      
      // Duplicates (should be rare). HEURISTICS here!
      if ( Abs( hit1.StartOnTarget( ) - hit2.StartOnTarget( ) ) < 5 ) {
	invalids[5] += 1;
	continue;
      }
      
      // A keeper.
      invalids[0] += 1;
      hit1.PrintParseable( best_out );
      hit2.PrintParseable( best_out );

      // Compute sep.
      const look_align &hit_rc = hit1.Rc1( ) ? hit1 : hit2;
      const look_align &hit_fw = hit1.Fw1( ) ? hit1 : hit2;
      int sep = hit_fw.StartOnTarget( ) - hit_rc.EndOnTarget( );
//      if ( CIRCULAR ) sep = Min( sep, sep + (int)hit_fw.TargetLength( ) );
      if ( CIRCULAR ) sep = ( sep + (int)hit_fw.TargetLength( ) ) % (int)hit_fw.TargetLength( ) ;
      seps_out << sep << "\n";
    }
    best_out.close( );
    seps_out.close( );

    // Generate table of discarded pairs.
    table.clear( );
    
    table.push_back( MkVec( String( " pairs in input" ),
			    ToStringAddCommas( BigSum( invalids ) ) ) );
    table.push_back( MkVec( String( " one or both end read unaligned" ),
			    ToStringAddCommas( invalids[1] ) ) );
    table.push_back( MkVec( String( " one or both end read multiply aligned" ),
			    ToStringAddCommas( invalids[2] ) ) );
    table.push_back( MkVec( String( " reads on different targets" ),
			    ToStringAddCommas( invalids[3] ) ) );
    table.push_back( MkVec( String( " wrong relative orientation" ),
			    ToStringAddCommas( invalids[4] ) ) );
    table.push_back( MkVec( String( " undetected duplicates" ),
			    ToStringAddCommas( invalids[5] ) ) );
    table.push_back( MkVec( String( " pairs left in output" ),
			    ToStringAddCommas( invalids[0] ) ) );
    
    log << "\n"
	<< "TABLE OF DISCARDED PAIRS\n"
	<< "\n";
    PrintTabular( log, table, 3, "rr" );
    log << endl;
    
  }

  // Build histogram of separations.
  histogram<int> histo;

  vec<int> bins;
  int step = 0;
  while( 1 ) {
    int bin = MIN_SEP + step * BIN_SIZE;
    if ( bin > MAX_SEP ) break;
    bins.push_back( bin );
    step++;
  }
  histo.AddBins( bins );
  
  int unders = 0;
  int overs = 0;
  int nval = 0;
  ifstream seps_in( seps_file.c_str( ) );
  while ( seps_in ) {
    int sep;
    seps_in >> sep;
    if ( !seps_in ) break;
    if ( sep < MIN_SEP ) unders++;
    else if ( sep > MAX_SEP ) overs++;
    else { nval++; histo.AddDatum( sep ); }
  }
  seps_in.close( );

  String str_nval = ToStringAddCommas( nval );
  String str_min = ToStringAddCommas( MIN_SEP );
  String str_max = ToStringAddCommas( MAX_SEP );
  String m1 = "Histogram of " + str_nval + " valid separations";
  String m2 = "Count of invalid separations (not in this plot):";
  String m3 = ToStringAddCommas( unders ) + " (because < " + str_min + ")";
  String m4 = ToStringAddCommas( overs ) + " (because > " + str_max + ")";
  
  ofstream txt_out( histo_txt_file.c_str( ) );
  txt_out << m1 << "\n"
	  << m2 << "\n"
	  << m3 << "\n"
	  << m4 << "\n"
	  << "\n";
  histo.PrintAsColumns( txt_out, true, true, false );
  txt_out << endl;
  txt_out.close( ) ;

  {
    using namespace ns_psplot;
    
    vec<ns_psplot::freetext> labels( 4 );
    labels[0] = freetext( m1, black, 16 );
    labels[1] = freetext( m2, black, 12 );
    labels[2] = freetext( m3, black, 12 );
    labels[3] = freetext( m4, black, 12 );

    ofstream eps_out( histo_eps_file.c_str( ) );
    histo.PrintAsEps( eps_out, labels, 0 );
    eps_out.close( );
  }
  
  // Done.
  String date = Date( );
  cout << date << ": ValidateFindIndexes done" << endl;
  log << date << ": ValidateFindIndexes done" << endl;
  log.close( );

}
