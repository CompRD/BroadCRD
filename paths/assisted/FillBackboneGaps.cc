///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <omp.h>   // needed by CInsertsDB
#include "MainTools.h"
#include "Basevector.h"
#include "PairsManager.h"
#include "SupersHandler.h"
#include "VecUtilities.h"
#include "feudal/IncrementalWriter.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"
#include "paths/KmerPathInterval.h"
#include "paths/ReadLoc.h"
#include "paths/assisted/CInsertsDB.h"
#include "paths/assisted/CloudWalker.h"
#include "paths/assisted/KmerPathPerfectOverlaps.h"
#include "util/RunCommand.h"
// MakeDepend: library OMP

/**
 * FillBackboneGaps
 *
 * Fill gaps in the given backbones, by walking a given global cloud
 * of reads / unibases (see UnifiedGlobalCloud). A smaller k-mer size
 * (SMALL_K) is used to merge partial closures, if found.  The main
 * args mirror those of RunAllPathsLG.
 */
int main( int argc, char *argv[] )
{
  RunTime( );
 
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String( SUBDIR );
  CommandArgument_String( ASSEMBLY_IN );
  CommandArgument_String( ASSEMBLY_OUT );

  // kmer size (big and small).
  CommandArgument_Int_OrDefault( K, 96 );
  CommandArgument_Int_OrDefault( SMALL_K, 24 );
  
  // Used to define valid gap sizes (with respect to input). Must be > 1.
  CommandArgument_Double_OrDefault( GAP_MULTIPLIER, 3.5 );

  // Used to define valid extension (and to limit search space).
  CommandArgument_Int_OrDefault( MIN_SEP, 50 );
  CommandArgument_Int_OrDefault( MAX_SEP, 15000 );
  CommandArgument_Int_OrDefault( MAX_KLEN, 5000 );

  // Max depth (for the recursive search).
  CommandArgument_Int_OrDefault( MAX_DEPTH, 512 );

  // Max dist (for secondary cloud, >=0. See CloudWalker.h for details).
  CommandArgument_Int_OrDefault( MAX_DIST, 0 );

  // Seed for randomizer.
  CommandArgument_Int_OrDefault( SEED, 666 );

  // Various core input files (use UnifiedGlobalCloud to build GLOBAL_CLOUD).
  CommandArgument_String_OrDefault( GLOBAL_CLOUD, "global_cloud" );
  CommandArgument_String_OrDefault( UNIPATHS, "all_reads" );
  CommandArgument_String_OrDefault( JUMPS, "jump_reads_ec" );

  // Log files and extra "saves" (some are optionals, for debugging).
  CommandArgument_Bool_OrDefault( VERBOSE_LOG, False );
  CommandArgument_Bool_OrDefault( DUMP_ALL_INSERTS, False );
  CommandArgument_Bool_OrDefault( DUMP_BRIEF_INSERTS, False );
  CommandArgument_Bool_OrDefault( SAVE_FILLED_ASSEMBLY, True );
  
  // Optionally, select a single gap for closure (by super id, contig pos).
  CommandArgument_Int_OrDefault( SELECT_SID, -1 );
  CommandArgument_Int_OrDefault( SELECT_SPOS, -1 );

  // Dump all full closures found (arg. to SelectClosures).
  CommandArgument_Bool_OrDefault( ALL_FULL, False );  

  // Thread control (for CInsertDB).
  CommandArgument_UnsignedInt_OrDefault( NUM_THREADS, 0 );

  EndCommandArguments;
  
  // Validate args.
  if ( GAP_MULTIPLIER < 1 ) {
    cout << "Fatal error: GAP_MULTIPLIER must be >= 1.\n" << endl;
    return 1;
  }

  // Thread control.
  NUM_THREADS = configNumThreads( NUM_THREADS );
  omp_set_num_threads( NUM_THREADS );

  // Dir and file names.
  String strK = ToString( K );
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  
  String assembly_in_head = sub_dir + "/" + ASSEMBLY_IN;
  String locs_head = assembly_in_head + ".contigs";
  String orig_ids_file = assembly_in_head + ".orig.ids";
  String locs_file = assembly_in_head + ".contigs.readlocs";
  String contigs_file = assembly_in_head + ".contigs.fastb";
  String supers_file = assembly_in_head + ".superb";
  
  String cloud_head = run_dir + "/" + GLOBAL_CLOUD;
  String to_cloud_map_file = cloud_head + ".to_gcid";
  String to_rc_cloud_paths_file = cloud_head + ".to_rc_paths.k" + strK;
  String cloud_paths_file = cloud_head + ".paths.k" + strK;
  String cloud_pathsdb_file = cloud_head + ".pathsdb.k" + strK;

  String unipaths_head = run_dir + "/" + UNIPATHS;
  String unipaths_file = unipaths_head + ".unipaths.k" + strK;
  String unipaths_rc_file = unipaths_head + ".unipaths_rc.k" + strK;
  String unipathsdb_file = unipaths_head + ".unipathsdb.k" + strK;
  String unibases_file = unipaths_head + ".unibases.k" + strK;

  String jumps_head = run_dir + "/" + JUMPS;
  String jumps_paths_file = jumps_head + ".paths.k" + strK;
  String jumps_pairs_file = jumps_head + ".pairs";
  
  String assembly_out_head = sub_dir + "/" + ASSEMBLY_OUT;
  String closures_head = assembly_out_head + ".closures";
  String main_log_file = assembly_out_head + ".FillBackboneGaps.log";
  String closures_fastb_file = closures_head + ".fastb";
  String closures_info_file = closures_head + ".info";
  String all_inserts_file = assembly_out_head + ".all_inserts.info";
  String brief_inserts_file = assembly_out_head + ".brief_inserts.info";

  vec<String> needed;
  needed.push_back( orig_ids_file );
  needed.push_back( locs_file );
  needed.push_back( contigs_file );
  needed.push_back( supers_file );
  needed.push_back( to_cloud_map_file );
  needed.push_back( to_rc_cloud_paths_file );
  needed.push_back( cloud_paths_file );
  needed.push_back( cloud_pathsdb_file );
  needed.push_back( unipaths_file );
  needed.push_back( unipaths_rc_file );
  needed.push_back( unipathsdb_file );
  needed.push_back( unibases_file );
  needed.push_back( jumps_paths_file );
  needed.push_back( jumps_pairs_file );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;

  // Check args.
  if ( SELECT_SID > -1 || SELECT_SPOS > -1 ) {
    if ( SELECT_SID == -1 || SELECT_SPOS == -1 ) {
      cout << "Fatal error, must set both SELECT_SID and SELECT_POS to fill\n"
	   << "a single gap. Leaving now.\n" << endl;
      return 1;
    }
  }

  // Log stream.
  ofstream log( main_log_file.c_str( ) );
  PrintCommandPretty ( log );
  cout << "Sending log to " << main_log_file << "\n" << endl;

  // Single gap selected.
  if ( SELECT_SID > -1 || SELECT_SPOS > -1 ) {
    log << "NB: filling only gap " << SELECT_SPOS
	<< " from super " << SELECT_SID
	<< "\n" << endl;
  }
  
  // Load.
  read_locs_on_disk locs_parser( locs_head, run_dir );

  log << Date( ) << ": loading cloud reads" << endl;
  vec<int> to_rc;
  vecKmerPath cloud_paths( cloud_paths_file );
  BinaryReader::readFile( to_rc_cloud_paths_file, &to_rc );
  READ( to_cloud_map_file, vec<int>, to_gcid );
  BREAD2( cloud_pathsdb_file, vec<tagged_rpint>, cloud_pathsdb );
  SCloud scloud( to_rc, cloud_paths, cloud_pathsdb );
  
  log << Date( ) << ": loading kbb" << endl;
  BREAD2( unipathsdb_file, vec<tagged_rpint>, unipathsdb );
  vecbvec unibases( unibases_file );
  vecKmerPath unipaths( unipaths_file );
  vecKmerPath unipaths_rc( unipaths_rc_file );
  KmerBaseBroker kbb( K, unipaths, unipaths_rc, unipathsdb, unibases );

  log << Date( ) << ": loading jump reads" << endl;
  vecKmerPath jpaths( jumps_paths_file );
  PairsManager jpairs( jumps_pairs_file );
  longlong njumps = jpaths.size( );
  
  log << Date( ) << ": loading supers" << endl;
  size_t n_contigs = MastervecFileObjectCount( contigs_file );
  shandler supers( n_contigs, supers_file );
  READ( orig_ids_file, vec<int>, orig_ids );
  
  log << Date( ) << ": loading orig ids" << endl;

  log << Date( ) << ": done loading\n" << endl;

  // Consts.
  const int infty = numeric_limits<int>::max( );

  // Initial CInsertDB (no super specified).
  int current_super_id = -1;
  CInsertsDB ins_db( MIN_SEP, MAX_SEP, 0, &supers, &locs_parser );
  
  // Log and out streams.
  ofstream all_inserts_out( DUMP_ALL_INSERTS
			    ? all_inserts_file.c_str( )
			    : "/dev/null" );
  
  ofstream brief_inserts_out( DUMP_BRIEF_INSERTS
			      ? brief_inserts_file.c_str( )
			      : "/dev/null" );
  
  ofstream cl_info_out( closures_info_file.c_str( ) );

  IncrementalWriter<BaseVec> cl_bases_out( closures_fastb_file.c_str( ) );
  
  // Loop over all supers.
  for (int super_id=0; super_id<supers.Size( ); super_id++) {
    if ( SELECT_SID > -1 && SELECT_SID != super_id ) continue;
    
    // Set super of CInsertDB.
    const superb &sup = supers[super_id];
    const int ngaps = sup.Ntigs( ) - 1;
    ins_db.SetSuper( &super_id );
    
    // Collect all inserts potentially filling gaps in the super.
    vec< pair<int,int> > cpairs;
    vec< vec<int> > iids;
    {
      vec< triple<int,int,int> > cpos2inserts = ins_db.CposToInserts( );

      if ( DUMP_ALL_INSERTS ) {
	vec< vec<String> > table;
	for (int jj=0; jj<cpos2inserts.isize( ); jj++)
	  table.push_back( ins_db.LineInfoAlt( cpos2inserts[jj].third ) );
	PrintTabular( all_inserts_out, table, 3, "rrrr" );
	all_inserts_out << "\n";
      }
      
      vec<int> emptyvec;
      for (int jj=0; jj<cpos2inserts.isize( ); jj++) {
	const triple<int,int,int> &trip = cpos2inserts[jj];
	if ( cpairs.size( ) < 1 ||
	     trip.first != cpairs.back( ).first ||
	     trip.second != cpairs.back( ).second ) {
	  cpairs.push_back( make_pair( trip.first, trip.second ) );
	  iids.push_back( emptyvec );
	}
	iids[iids.size( )-1].push_back( cpos2inserts[jj].third );
      }

      SortSync( cpairs, iids );
    }
    
    // Cloud walk all of these.
    for (int gap_id=0; gap_id<ngaps; gap_id++) {
      if ( SELECT_SPOS > -1 && SELECT_SPOS != gap_id ) continue;
      
      // Log event.
      int gap = sup.Gap( gap_id );
      log << "Walking gap between s"
	  << super_id << "." << gap_id << "/" << sup.Ntigs( )
	  << " = c" << sup.Tig( gap_id ) << " and s"
	  << super_id << "." << gap_id + 1 << "/" << sup.Ntigs( )
	  << " = c" << sup.Tig( gap_id + 1 )
	  << " (gap size: " << gap << " +/- "
	  <<  sup.Dev( gap_id ) << ")\n";
      
      // Pos of flanking contigs.
      pair<int,int> target_pair = make_pair( gap_id, gap_id + 1 );
      
      // Build list of candidates for walking.
      vec<int> iids_sel;
      vec< pair<int,int> >::iterator it;
      it = lower_bound( cpairs.begin( ), cpairs.end( ), target_pair );
      
      // All links between contigs at target_pair.
      while ( it != cpairs.end( ) ) {
	if ( *it != target_pair ) break;
	const vec<int> &all_ids = iids[ distance( cpairs.begin( ), it ) ];
	iids_sel.reserve( all_ids.size( ) );
	for (int jj=0; jj<all_ids.isize( ); jj++) {
	  if ( gap_id == 0 || gap_id == ngaps - 1 ) {
	    if ( ins_db.IsSeparatedType( all_ids[jj] ) )
	      iids_sel.push_back( all_ids[jj] );
	  }
	  else {
	    if ( ins_db.IsValidType( all_ids[jj] ) )
	      iids_sel.push_back( all_ids[jj] );
	  }
	}
	it++;
      }

      // Report links count.
      if ( DUMP_BRIEF_INSERTS ) {
	brief_inserts_out << "s" << super_id << "\t"
			  << "gap " << target_pair.first
			  << " to " << target_pair.second
			  << "\t" << iids_sel.size( )
			  << " candidates" << "\n";
	if ( gap_id == ngaps - 1 ) brief_inserts_out << "\n";
      }
      
      // Indexes: contig_id -> unibase_id -> global_cloud_id.
      int u1 = orig_ids[ sup.Tig( target_pair.first ) ];
      int u2 = orig_ids[ sup.Tig( target_pair.second ) ];
      bool u1rc = ( u1 < 0 );
      bool u2rc = ( u2 < 0 );
      if ( u1rc ) u1 = - 1 - u1;
      if ( u2rc ) u2 = - 1 - u2;
      int gcid1 = u1rc ? to_rc[ to_gcid[u1] ] : to_gcid[u1];
      int gcid2 = u2rc ? to_rc[ to_gcid[u2] ] : to_gcid[u2];
      
      int u1_len = unibases[u1].size( );
      int u2_len = unibases[u2].size( );
      int pos1 = target_pair.first;
      int pos2 = target_pair.second;
      ForceAssertEq( u1_len, sup.Len( pos1 ) );
      ForceAssertEq( u2_len, sup.Len( pos2 ) );

      int gclen = cloud_paths[gcid1].TotalLength( ) + K - 1;
      ForceAssertEq( gclen, sup.Len( pos1 ) );
      
      // Find closures.
      CClosures closures;
      FillGap( SMALL_K, K, MAX_KLEN, MAX_DEPTH, MAX_DIST, VERBOSE_LOG, gap_id,
	       gcid1, gcid2, scloud, ins_db, jpaths, kbb, closures, &log );
      
      // Chomp K-1 bases from left and/or right side of closures.
      closures.ChompAdjacencies( K );
      
      // Sort closures.
      closures.GapBasedSort( gap );

      // Closures found, select winners.
      SelectClosures( K, gap, GAP_MULTIPLIER, ALL_FULL, closures );

      // No closures found.
      if ( closures.Size( ) < 1 ) continue;
      
      // Closures found: dump closures and some info.
      for (int ii=0; ii<(int)closures.Size( ); ii++) {
	cl_bases_out.add( closures.Bases( ii ) );
	closures.PrintInfo( ii, super_id, target_pair.first, cl_info_out );
      }
      
    } // loop over all gaps in super
    
  } // loop over all supers
  
  // Close streams.
  if ( DUMP_ALL_INSERTS ) all_inserts_out.close( );
  if ( DUMP_BRIEF_INSERTS ) brief_inserts_out.close( );
  cl_info_out.close( );
  cl_bases_out.close( );

  // Save assembly.
  if ( SAVE_FILLED_ASSEMBLY ) {
    log << Date( ) << ": saving" << endl;
    SaveFilledAssembly( closures_head, assembly_in_head, assembly_out_head );
  }
  else
    log << Date( ) << ": NOTE: filled assembly not saved" << endl;
  
  // Done.
  String str_date = Date( );
  cout << str_date << ": done" << endl;
  log << str_date << ": done" << endl;

}

