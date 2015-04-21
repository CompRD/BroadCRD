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
#include "PairsManager.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"
#include "lookup/QueryLookupTableCore.h"
#include "math/Functions.h"
#include "random/Shuffle.h"

/**
 * EvalLibOnUnibases
 *
 * Select <SAMPLE_SIZE> pairs from a given library, and align the
 * reads to the specified <UNIBASES> (after stripping the last <K>-1
 * bases from the unibases). The alignments are minimally filtered.
 * Output is sent to the main log file in <READS>.<OUTBASE>.
 * 
 * READS: it loads <READS>.{fastb,pairs}
 * UNIBASES: it loads <UNIBASES>.unibases.k<K>
 * OUTBASE: output is saved in <READS>.<OUTBASE>
 * MIN_LENGTH: min length to align (after stripping last k-1 bases)
 * SAMPLE_SIZE: how many pairs to align (randomly selected)
 * MAX_MUTATION_RATE: discard aligns with excessive mutation rate
 * MAX_INDELS: discard aligns with too many indels
 * SEED: used to seed the randomizer
 * FORCE: do not use cached data, even if found
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_Int( K );
  CommandArgument_String( READS );
  CommandArgument_String( UNIBASES );
  CommandArgument_String_OrDefault( OUTBASE, "EvalLib" );
  CommandArgument_Int_OrDefault( MIN_LENGTH, 2000 );
  CommandArgument_Int_OrDefault( SAMPLE_SIZE, 100000 );
  CommandArgument_Double_OrDefault( MAX_MUTATION_RATE, 0.15 );
  CommandArgument_Int_OrDefault( MAX_INDELS, 2 );
  CommandArgument_Int_OrDefault( SEED, 666 );
  CommandArgument_Bool_OrDefault( FORCE, True );
  EndCommandArguments;
  
  // Dir and file names.
  String strK = ToString( K );
  
  String out_dir = READS + "." + OUTBASE;

  String reads_fastb_file = READS + ".fastb";
  String reads_pairs_file = READS + ".pairs";
  String unibases_file = UNIBASES + ".unibases.k" + strK;
  String log_file = out_dir + "/main.log";
  String target_head = out_dir + "/target";
  String lookup_file = target_head + ".lookup";
  String target_file = target_head + ".fastb";
  String tids_file = target_head + ".ids";
  String query_file = out_dir + "/query.fastb";
  String qids_file = out_dir + "/query.ids";
  String all_hits_file = out_dir + "/all_aligns.qlt";
  String filtered_hits_file = out_dir + "/filtered_aligns.qlt";
  
  Mkpath( out_dir );
  ofstream log( log_file.c_str( ) );
  PrintCommandPretty( log );
  cout << "See " << log_file << " for details\n" << endl;
  
  // Load.
  log << Date( ) << ": loading pairing info" << endl;
  PairsManager pairs( reads_pairs_file );
  uint64_t n_pairs = pairs.nPairs( );

  // Digest unibases, and save target file. Generate lookup.
  if ( FORCE || ! IsRegularFile( lookup_file ) ) {
    log << Date( ) << ": loading unibases" << endl;
    vecbvec unibases_orig( unibases_file );

    log << Date( ) << ": digesting unibases" << endl;
    vec<int> tids;
    vecbvec unibases;
    tids.reserve( unibases_orig.size( ) );
    unibases.reserve( unibases_orig.size( ) );
    for (size_t id=0; id<unibases_orig.size( ); id++) {
      if ( unibases_orig[id].size( ) - (K-1) < size_t(MIN_LENGTH) ) continue;
      int len = unibases_orig[id].size( ) - ( K - 1 );
      tids.push_back( id );
      unibases.push_back( bvec( unibases_orig[id], 0, len ) );
    }
    
    log << "\n"
	<< "Unibases in input: " << unibases_orig.size( ) << "\n"
	<< "unibases selected: " << unibases.size( ) << "\n"
	<< "(after removal of unibases shorter than " << MIN_LENGTH << " bp)\n"
	<< endl;
    
    log << Date( ) << ": saving unibases" << endl;
    WRITE( tids_file,tids );
    unibases.WriteAll( target_file );

    log << Date( ) << ": generating lookup table" << endl;
    String comm
      = "MakeLookupTable SOURCE=" + target_file
      + " OUT_HEAD=" + target_head
      + " LOOKUP_ONLY=True";
    SystemSucceedQuiet( comm );
  }    

  // Randomly select pairs, and save query file.
  if ( FORCE || ! IsRegularFile( query_file ) ) {
    vec<uint64_t> shuffled;
    Shuffle64( n_pairs, shuffled, uint64_t( SEED ) );
    
    vec<int> qids;   // WARNING! Must use int for SparseRead to work!
    uint64_t n_sel = Min( n_pairs, uint64_t( SAMPLE_SIZE ) );
    qids.reserve( 2 * n_sel );
    for (uint64_t ii=0; ii<n_sel; ii++) {
      qids.push_back( (int)pairs.ID1( shuffled[ii] ) );
      qids.push_back( (int)pairs.ID2( shuffled[ii] ) );
    }
    sort( qids.begin( ), qids.end( ) );
    
    log << "\n"
	<< "Selected " << qids.size( )
	<< " reads for alignment\n"
	<< endl;

    log << Date( ) << ": loading selected reads" << endl;
    vecbvec query;
    query.SparseRead( reads_fastb_file, qids, 0 );
    
    vecbvec qselect;
    qselect.reserve( qids.size( ) );
    for (size_t ii=0; ii<qids.size( ); ii++)
      qselect.push_back( query[ qids[ii] ] );
    
    log << Date( ) << ": saving selected reads" << endl;
    WRITE( qids_file, qids );
    qselect.WriteAll( query_file );
  }
  
  // Align reads.
  String qlt_args = "K=12 MM=12 MF=5000 SH=True MC=0.15";
  if ( FORCE || ! IsRegularFile( all_hits_file ) ) {
    log << Date( ) << ": aligning reads" << endl;
    String qlt_comm
      = "QueryLookupTable"
      + String( " AI=True" )
      + String( " PARSEABLE=True" )	    
      + String( " PRINT_MEMORY_USAGE=True" )
      + String( " LIST_UNPLACED_BY_PASS=False" )
      + String( " " + qlt_args )
      + String( " L=" + lookup_file )
      + String( " SEQS=" + query_file )
      + String( " TMP_DIR=" + out_dir )
      + String( " >& " + all_hits_file );
    SystemSucceed( qlt_comm );
  }

  // Filter aligns or load filtered aligns.
  vec<look_align_plus> hits;
  if ( FORCE || ! IsRegularFile( filtered_hits_file ) ) {
    log << Date( ) << ": loading aligns" << endl;
    vec<look_align_plus> all_hits;
    LoadLookAlignPlus( all_hits_file, all_hits );

    log << Date( ) << ": filtering aligns" << endl;
    hits.reserve( all_hits.size( ) );
    for (size_t ii=0; ii<all_hits.size( ); ii++) {
      const look_align_plus &hit = all_hits[ii];

      // Read is not fully embedded.
      if ( hit.pos1( ) > 0 || hit.Pos1( ) < (int)hit.QueryLength( ) ) continue;
      
      // Too many mismatches.
      if ( hit.MutationRate( ) > MAX_MUTATION_RATE ) continue;

      // Too many indels.
      if ( hit.indels > MAX_INDELS ) continue;

      // Ok.
      hits.push_back( hit );
    }
    
    log << Date( ) << ": saving filtered aligns" << endl;
    ofstream out( filtered_hits_file.c_str( ) );
    for (size_t ii=0; ii<hits.size( ); ii++)
      hits[ii].WriteParseable( out );
    out.close( );
  }
  else {
    log << Date( ) << ": loading filtered aligns" << endl;
    LoadLookAlignPlus( filtered_hits_file, hits );
  }

  // Early exit (no aligns found).
  if ( hits.size( ) < 1 ) {
    log << "\nNo aligns found.\n\n" << Date( ) << ": done" << endl;
    log.close( );
    return 0;
  }

  // Final stats.
  log << Date( ) << ": generating summary statistics" << endl;
  vecbvec target( target_file );
  longlong total_tlen = 0;
  for (size_t ii=0; ii<target.size( ); ii++)
    total_tlen += target[ii].size( );

  size_t n_tested = MastervecFileObjectCount( query_file );
  vec<bool> aligned( n_tested, 0 );
  for (size_t ii = 0; ii < hits.size(); ii++)
    aligned[hits[ii].query_id] = True;
  size_t n_aligned = 0;
  for (size_t ii=0; ii<n_tested; ii++)
    if ( aligned[ii] ) n_aligned++;

  double ratio = double( n_aligned ) / double( n_tested );

  longlong read_len = hits[0].query_length;
  double cov1 = double( read_len * n_aligned ) / double( total_tlen );

  size_t n_total = MastervecFileObjectCount( reads_fastb_file );
  double cov2 = cov1 * double( n_total ) / double(  n_tested );

  double ratio2 = double( n_tested ) / double( n_total );

  // Print info (for now no pairing info).
  log << "\n"
      << "FINAL STATISTICS\n"
      << "\n"
      << "qlt arguments (heuristics): " << qlt_args << "\n"
      << "reads in this library: " << n_total << "\n"
      << "reads tested for alignment: " << n_tested << "\n"
      << "fraction tested: " << ToString( 100.0 * ratio2, 1 ) << "\n"
      << "reads aligned: " << n_aligned << "\n"
      << "fraction aligned (wrt tested): " << ToString(100.0* ratio,2) << "%\n"
      << "coverage of aligned reads: " << ToString( cov1, 1 ) << "x\n"
      << "estimated coverage of full set: " << ToString( cov2, 1 ) << "x\n"
      << "\n"
      << "NB: coverage is based on an estimated genome size of "
      << total_tlen << " bp,\n"
      << "as computed by summing up the lengths of unibases longer than"
      << " " << MIN_LENGTH << " bp\n"
      << "(as specified by the argument MIN_LENGTH).\n"
      << "\n";

  // Done.
  log << Date( ) << ": done" << endl;
  log.close( );

}
