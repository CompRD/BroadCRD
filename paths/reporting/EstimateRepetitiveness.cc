///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "VecUtilities.h"
#include "math/NStatsTools.h"
#include "paths/CommonPatherCore.h"
#include "util/RunCommand.h"
#include "feudal/BinaryStream.h"

/**
 * EstimateRepetitiveness
 *
 * Estimate repetitiveness on a given fastb by reporting N-stats on
 * unique vs repetitive kmers stretches. It only needs a fastb in
 * input.
 *
 * NUM_THREADS: used by CommonPather
 * FORCE: do not use cached data (and clean up tmp files at the end)
 * VERBOSE: verbose log flag
 */
int main( int argc, char *argv[] )
{

  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( FASTB );
  CommandArgument_String_OrDefault( TMP, FASTB + ".EstimateRepetitiveness" );
  CommandArgument_Int_OrDefault( K, 100 );
  CommandArgument_Int_OrDefault( NUM_THREADS, 8 );
  CommandArgument_Bool_OrDefault( FORCE, True );
  CommandArgument_Bool_OrDefault( VERBOSE, False );
  EndCommandArguments;
  
  // Dir and file names.
  String strK = "k" + ToString( K );
  
  string contigs_head = TMP + "/contigs";
  String fastb_file = contigs_head + ".fastb";
  String paths_file = contigs_head + ".paths";
  String pathsdb_file = contigs_head + ".pathsdb." + strK;
  String pathsdb_big_file = contigs_head + ".pathsdb_big." + strK;

  Mkpath( TMP );
  cout << "All lengths in kmers!\n" << endl;
  
  // Run CommonPather.
  if ( FORCE || ! IsRegularFile( fastb_file ) )
    Cp( FASTB, fastb_file );
  
  if ( FORCE || ! IsRegularFile( paths_file + "." + strK ) ) {
    vec<String> fastb_in( 1, fastb_file );
    vec<String> paths_out( 1, paths_file );
    CommonPather( K, fastb_in, paths_out, NUM_THREADS );
  }

  // Run MakeRcDb (logging within).
  bool db_exists = IsRegularFile( pathsdb_file );
  bool db_big_exists = IsRegularFile( pathsdb_big_file );
  if ( FORCE || ! ( db_exists || db_big_exists ) ) {
    String theCommand
      = "MakeRcDb PRE= DATA= RUN= K=" + ToString( K )
      + " READS=" + contigs_head;
    RunCommand( theCommand );
  }
  
  // Load paths, and pathsdb.
  vecKmerPath paths( paths_file + "." + strK );

  bool big_db = IsRegularFile( pathsdb_big_file );
  vec<tagged_rpint> pathsdb;
  vec<big_tagged_rpint> pathsdb_big;
  if ( big_db ) BinaryReader::readFile( pathsdb_big_file, &pathsdb_big );
  else BinaryReader::readFile( pathsdb_file, &pathsdb );
  
  // No paths found (probably an error).
  if ( paths.size( ) < 1 ) {
    cout << "No paths found. Exit\n" << endl;
    return 0;
  }

  // Lengths of cn1 and cn>1 intervals;
  vec<int> cn1_lens;
  vec<int> cn2_lens;   // all cn>1 in here
  size_t n_lens = 0;
  for (size_t ii=0; ii<paths.size( ); ii++)
    n_lens += paths[ii].NSegments( );
  cn1_lens.reserve( n_lens );
  cn1_lens.reserve( n_lens );

  // Lengths of paths.
  vec<int> path_lens;
  path_lens.reserve( paths.size( ) );

  // Loop over all paths.
  for (size_t path_id=0; path_id<paths.size( ); path_id++) {
    const KmerPath &kpath = paths[path_id];
    path_lens.push_back( kpath.TotalLength( ) );
    
    vec< pair<int,int> > len2cn;   // length of interval, cn of interval
    
    // Loop over all segments.
    for (int segment_id=0; segment_id<kpath.NSegments( ); segment_id++) {
      const KmerPathInterval &segment = kpath.Segment( segment_id );
      const int segment_len = kpath.Length( segment_id );
      
      // Loop over all kmer ids in the segment.
      for (longlong ii=segment.Start( ); ii<=segment.Stop( ); ii++) {
	vec<longlong> locs;      
	if ( big_db ) Contains( pathsdb_big, ii, locs, false, 3 );
	else Contains( pathsdb, ii, locs, false, 3 );
	ForceAssertGt( locs.isize( ), 0 );
	int cn = locs.size( ) == 1 ? 1 : 2;
	
	if ( len2cn.size( ) < 1 )
	  len2cn.push_back( make_pair( 1, cn ) );
	else {
	  bool append = ( len2cn.back( ).second == cn );
	  if ( append ) len2cn[len2cn.size( )-1].first += 1;
	  else len2cn.push_back( make_pair( 1, cn ) );
	}
      }
    }    

    // Update cn1_lens, cn2_lens.
    for (int ii=0; ii<len2cn.isize( ); ii++) {
      if ( len2cn[ii].second == 1 ) cn1_lens.push_back( len2cn[ii].first );
      else cn2_lens.push_back( len2cn[ii].first );
    }    

    // Print verbose info.
    if ( ! VERBOSE ) continue;

    int begin = 0;
    vec< vec<String> > table;
    table.push_back( MkVec( String( "contig" ),
			    String( "kmer_length" ),
			    String( "begin" ),
			    String( "end" ),
			    String( "" ) ) );
    for (size_t ii=0; ii<len2cn.size( ); ii++) {
      int len = len2cn[ii].first;
      String tag = len2cn[ii].second == 1 ? "" : "REPEAT";
      table.push_back( MkVec( ToString( path_id ),
			      ToString( len ),
			      ToString( begin ),
			      ToString( begin + len ),
			      tag ) );
      begin += len;
    }
    PrintTabular( cout, table, 3, "rrrrl" );
    cout << endl;

  }

  // Print stats (vectors will be sorted).
  vec< vec<String> > table;
  longlong tot_pathlen = BigSum( path_lens );
  longlong tot_unique = BigSum( cn1_lens );
  longlong tot_rep = BigSum( cn2_lens );
  double pc_unique = 100.0 * SafeQuotient( tot_unique, tot_pathlen );
  double pc_rep = 100.0 * SafeQuotient( tot_rep, tot_pathlen );
  table.push_back( MkVec( String( "total path length" ),
			  ToStringAddCommas( tot_pathlen ),
			  String( "" ) ) );
  table.push_back( MkVec( String( "total unique length" ),
			  ToStringAddCommas( tot_unique ),
			  ToString( pc_unique, 2 ) + "%" ) );
  table.push_back( MkVec( String( "total repetitive length" ),
			  ToStringAddCommas( tot_rep ),
			  ToString( pc_rep, 2 ) + "%" ) );
  PrintTabular( cout, table, 3, "lrr" );
  cout << endl;
  
  PrintBasicNStats( "paths", path_lens, cout );
  cout << endl;

  PrintBasicNStats( "unique", cn1_lens, cout );
  cout << endl;

  PrintBasicNStats( "repetitive", cn2_lens, cout );
  cout << endl;
  
  // Clean up.
  if ( FORCE ) {
    String theCommand = "rm -rf " + TMP;
    RunCommandWithLog( theCommand, "/dev/null" );
  }

  // Done.
  cout << Date( ) << ": done" << endl;

}

