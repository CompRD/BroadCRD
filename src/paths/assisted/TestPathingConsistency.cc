///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "PairsManager.h"
#include "Superb.h"
#include "VecUtilities.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"
#include "paths/KmerPathInterval.h"
#include "util/RunCommand.h"

void TestReads( const int &K,
		const vecKmerPath &reads_paths,
		const vecbvec &reads_bases,
		const vec<tagged_rpint> &unipathsdb,
		const KmerBaseBroker &kbb,
		const String &descriptor,
		ostream &log );

/**
 * TestPathingConsistency
 *
 * Kmer consistency tests for the path based assisted assembly code.
 * Make sure reads (frags and/or jumps) live in the same kmer space as
 * a given set of unibases.
 *
 * HEAD: head name of assembly
 * CONTIGS: optional extra ".contigs" before .fastb (may be "")
 * UNIPATHS: head name relative to RUN
 * FRAGS: head name of frag reads, relative to RUN
 * JUMPS: head name of jump reads, relative to RUN
 */
int main( int argc, char *argv[] )
{
  RunTime( );
 
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String( SUBDIR );
  CommandArgument_String( HEAD );
  CommandArgument_Int_OrDefault( K, 96 );
  CommandArgument_String_OrDefault( CONTIGS, ".contigs" );
  CommandArgument_String_OrDefault( UNIPATHS, "all_reads" );
  CommandArgument_String_OrDefault( FRAGS, "filled_reads_filt");
  CommandArgument_String_OrDefault( JUMPS, "jump_reads_ec");
  EndCommandArguments;
  
  // Dir and file names.
  String strK = ToString( K );
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  
  String assembly_head = sub_dir + "/" + HEAD + ".k" + strK;
  String to_orig_map_F = assembly_head + ".orig.ids";
  String contigs_F = assembly_head + CONTIGS + ".fastb";

  String unipaths_head = run_dir + "/" + UNIPATHS;
  String frags_head = run_dir + "/" + FRAGS;
  String jumps_head = run_dir + "/" + JUMPS;
  
  String unipaths_F = unipaths_head + ".unipaths.k" + strK;
  String unipaths_rc_F = unipaths_head + ".unipaths_rc.k" + strK;
  String unipathsdb_F = unipaths_head + ".unipathsdb.k" + strK;
  String unibases_F = unipaths_head + ".unibases.k" + strK;
  
  String frags_paths_F = frags_head + ".paths.k" + strK;
  String frags_bases_F = frags_head + ".fastb";
  bool use_frags = ( FRAGS != "" );

  String jumps_paths_F = jumps_head + ".paths.k" + strK;
  String jumps_bases_F = jumps_head + ".fastb";
  bool use_jumps = ( JUMPS != "" );

  if ( ! ( use_frags || use_jumps ) ) {
    cout << "No frags or jumps specified. Leaving now.\n" << endl;
    return 1;
  }

  vec<String> needed;
  needed.push_back( to_orig_map_F );
  needed.push_back( contigs_F );
  needed.push_back( unipaths_F );
  needed.push_back( unipathsdb_F );
  needed.push_back( unipaths_rc_F );
  if ( use_frags ) needed.push_back( frags_paths_F );
  if ( use_frags ) needed.push_back( frags_bases_F );
  if ( use_jumps ) needed.push_back( jumps_paths_F );
  if ( use_jumps ) needed.push_back( jumps_bases_F );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;

  // Load.
  cout << Date( ) << ": loading paths and unipaths" << endl;
  vecKmerPath unipaths( unipaths_F );
  vecKmerPath frags_paths( frags_paths_F );
  vecKmerPath jumps_paths( jumps_paths_F );

  cout << Date( ) << ": loading reads" << endl;
  vecbvec frags_bases;
  vecbvec jumps_bases;
  if ( use_frags ) frags_bases.ReadAll ( frags_bases_F );
  if ( use_jumps ) jumps_bases.ReadAll( jumps_bases_F );
  
  cout << Date( ) << ": loading the unipaths' kbb" << endl;
  BREAD2( unipathsdb_F, vec<tagged_rpint>, unipathsdb );
  vecbvec unibases( unibases_F );
  vecKmerPath unipaths_rc( unipaths_rc_F );
  KmerBaseBroker kbb( K, unipaths, unipaths_rc, unipathsdb, unibases );

  cout << Date( ) << ": loading contigs" << endl;
  vecbvec contigs( contigs_F );

  // Test reads.
  if ( use_frags )
    TestReads( K, frags_paths, frags_bases, unipathsdb, kbb, "frag", cout );

  if ( use_jumps )
    TestReads( K, jumps_paths, jumps_bases, unipathsdb, kbb, "jump", cout );

  // Done.
  cout << Date( ) << ": done" << endl;

}

/**
 * TestReads
 *
 * Make sure kmers in reads are in unipaths, and check bases from the
 * fastb of the reads match with the bases of the kmers from the
 * KmerBaseBroker.
 */
void TestReads( const int &K,
		const vecKmerPath &reads_paths,
		const vecbvec &reads_bases,
		const vec<tagged_rpint> &unipathsdb,
		const KmerBaseBroker &kbb,
		const String &descriptor,
		ostream &log )
{
  longlong dotter = 1000000;
  longlong n_tested = 0;
  longlong n_failedA = 0;
  longlong n_failedB = 0;

  log << Date( ) << ": parsing "
      << ToStringAddCommas( reads_paths.size( ) ) << " "
      << descriptor << " reads (. = "
      << ToStringAddCommas( dotter ) << " reads)" << endl;
  for (size_t read_id=0; read_id<reads_paths.size( ); read_id++) {
    if ( read_id % dotter == 0 ) Dot( log, read_id / dotter );
    const bvec &bases = reads_bases[read_id];
    const KmerPath &path = reads_paths[read_id];

    int pos = 0;
    for (int segment_id=0; segment_id<path.NSegments( ); segment_id++) {
      const KmerPathInterval &seg = path.Segment( segment_id );
      for (longlong kmer=seg.Start( ); kmer<=seg.Stop( ); kmer++) {
	n_tested++;
	pos++;
	vec<longlong> locs;
	Contains( unipathsdb, kmer, locs, false, 1 );

	// Look for kmer in unipathsdb.
	if ( locs.size( ) < 1 ) {
	  n_failedA++;
	  continue;
	}

	// Sanity check bases.
	const bvec &b1 = kbb.Bases( kmer );
	bvec b2( bases, pos-1, K );
	if ( b1 != b2 )
	  n_failedB++;
      }
    }
    
  }
  log << endl;

  log << "\n"
      << "OVERALL STATS FOR " << descriptor << " READS\n"
      << "\n"
      << "total:              " << ToStringAddCommas( n_tested ) << "\n"
      << "not in unipathsdb:  " << ToStringAddCommas( n_failedA ) << "\n"
      << "bases do not match: " << ToStringAddCommas( n_failedB ) << "\n"
      << endl;
}

