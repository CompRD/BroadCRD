/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// MakeDepend: dependency MakeRcDb
#include "CoreTools.h"
#include "MainTools.h"
#include "Superb.h"
#include "VecUtilities.h"
#include "paths/CommonPatherCore.h"
#include "paths/KmerPath.h"
#include "paths/PdfEntry.h"
#include "paths/UnibaseUtils.h"
#include "util/ReadTracker.h"
#include "util/RunCommand.h"
#include "feudal/BinaryStream.h"

void SingleGapExaminer ( const String &out_dir,
			 const String &contigs_paths_file,
			 const digraph &adjacency,
			 const VecPdfEntryVec &copynum,
			 const vecKmerPath &unipaths_paths,
			 const vec<tagged_rpint> &unipathsdb,
			 const vec<superb> &supers,
			 const int CONTIG_ID,
			 const int DEPTH,
			 const int MAX_CN,
			 const bool EPS );

/**
 * GapsExaminer
 *
 * Interactive tool to look at a gap between a contig and its
 * successor in a scaffold. It tries to "walk" over the gap using the
 * unipath graph. Output is saved as:
 *   ../<SUBDIR>/<OUTDIR>/gap_c<CONTIG_ID>.log
 *   ../<SUBDIR>/<OUTDIR>/gap_c<CONTIG_ID>.{dot,eps}   (optional)
 * where CONTIG_ID is provided interactively.
 *
 * Warning: the first time the code is run, it will generate (and
 * cache) several files (this is a time consuming part).
 *
 * OUT_DIR: relative to SUBDIR
 * FORCE: do not use cached files
 * EPS: generate a dot file (save both .dot and .eps)
 * DEPTH: max depth to walk on the graph
 * MAX_CN: stop walking at repetitive unipaths
 * NUM_THREADS: use all available processors if 0
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_Int( K );
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String_OrDefault( READS, "all_reads" );
  CommandArgument_String_OrDefault( ASSEMBLY, "final" );
  CommandArgument_String_OrDefault( OUT_DIR, ASSEMBLY + ".GapsExaminer" );
  CommandArgument_Bool_OrDefault( FORCE, False );
  CommandArgument_Bool_OrDefault( EPS, True );
  CommandArgument_Int_OrDefault( DEPTH, 12 );
  CommandArgument_Int_OrDefault( MAX_CN, 24 );
  CommandArgument_Int_OrDefault( NUM_THREADS, 0 );
  EndCommandArguments;

  // Dir and file names.
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  String out_dir = sub_dir + "/" + OUT_DIR;

  String strK = ToString( K );
  String reads_base = run_dir + "/" + READS;
  String unipaths_bases_file =  reads_base + ".unibases.k" + strK;
  String unipaths_cn_file = reads_base + ".unipaths.predicted_count.k" + strK;
  String contigs_bases_file = sub_dir + "/" + ASSEMBLY + ".contigs.fastb";
  String supers_file = sub_dir + "/" + ASSEMBLY + ".superb";
  String unipaths_adjacency_file = out_dir + "/unipaths.adjacency.k" + strK;
  String unipaths_paths_file = out_dir + "/unipaths.paths"; 
  String contigs_paths_file = out_dir + "/contigs.paths"; 
  String full_unipaths_paths_file = unipaths_paths_file + ".k" + strK;
  String full_unipaths_pathsdb_file = unipaths_paths_file + "db.k" + strK;
  String full_contigs_paths_file = out_dir + "/contigs.paths.k" + strK; 
  
  Mkpath( out_dir );
  
  // Thread controls.
  NUM_THREADS = configNumThreads( NUM_THREADS );
  
  // Run CommonPather on unipath graph plus assembly contigs.
  if ( FORCE || ! IsRegularFile( full_unipaths_paths_file ) ) {
    vec<String> bases_in = MkVec( unipaths_bases_file, contigs_bases_file );
    vec<String> paths_out = MkVec( unipaths_paths_file, contigs_paths_file );
    CommonPather( K, bases_in, paths_out, NUM_THREADS );
  }
  
  // Run MakeRcDb.
  if ( FORCE || ! IsRegularFile( full_unipaths_pathsdb_file ) ) {
    cout << Date( ) << ": running MakeRcDb" << endl;
    String theCommand
      = "MakeRcDb K=" + ToString( K )
      + " PRE= DATA= RUN= READS=" + out_dir + "/unipaths";
    RunCommandWithLog( theCommand, String( "/dev/null" ) );
  }
  
  // Load adjacency graph (build it if needed).
  digraph adjacency;
  if ( FORCE || ! IsRegularFile( unipaths_adjacency_file ) ) {
    cout << Date( ) << ": loading unibases" << endl;
    vecbvec unibases( unipaths_bases_file );
    
    cout << Date( ) << ": building adjacency graph" << endl;
    BuildUnibaseAdjacencyGraph( unibases, adjacency, K );

    cout << Date( ) << ": saving adjacency graph" << endl;
    BinaryWriter::writeFile( unipaths_adjacency_file, adjacency );
  }
  else {
    cout << Date( ) << ": loading adjacency graph" << endl;
    BinaryReader::readFile( unipaths_adjacency_file, &adjacency );
  }

  // Load predicted copy number.
  cout << Date( ) << ": loading predicted copy numbers" << endl;
  VecPdfEntryVec copynum( unipaths_cn_file );

  // Load unipaths.
  cout << Date( ) << ": loading unipaths.paths" << endl;
  vecKmerPath unipaths_paths( full_unipaths_paths_file );

  cout << Date( ) << ": loading unipaths.pathsdb" << endl;
  BREAD2( full_unipaths_pathsdb_file, vec<tagged_rpint>, unipathsdb );
  
  // Load supers, and find gaps in super.
  cout << Date( ) << ": loading superbs" << endl;
  vec<superb> supers;
  ReadSuperbs( supers_file, supers );
  
  // Interactively call SingleGapExaminer.
  while ( 1 ) {
    cout << "> " << flush;
    
    String cmnd;
    getline( cin, cmnd );
    if ( ! cin || cmnd == "q" ) {
      cout << endl;
      break;
    }

    if ( ! cmnd.IsInt( ) ) {
      cout << cmnd << " is not a CONTIG_ID." << endl;
      continue;
    }
    int CONTIG_ID = cmnd.Int( );
    
    SingleGapExaminer( out_dir, full_contigs_paths_file, adjacency,
		       copynum, unipaths_paths, unipathsdb, supers,
		       CONTIG_ID, DEPTH, MAX_CN, EPS );
  }
  
  // Done.
  cout << Date( ) << ": done" << endl;

}

/**
 * SingleGapExaminer
 */
void SingleGapExaminer( const String &out_dir,
			const String &contigs_paths_file,
			const digraph &adjacency,
			const VecPdfEntryVec &copynum,
			const vecKmerPath &unipaths_paths,
			const vec<tagged_rpint> &unipathsdb,
			const vec<superb> &supers,
			const int CONTIG_ID,
			const int DEPTH,
			const int MAX_CN,
			const bool EPS )
{
  // File names (and open log stream).
  String dot_file = out_dir + "/gap_c" + ToString( CONTIG_ID ) + ".dot";
  String eps_file = out_dir + "/gap_c" + ToString( CONTIG_ID ) + ".eps";
  String log_file = out_dir + "/gap_c" + ToString( CONTIG_ID ) + ".log";
  ofstream log( log_file.c_str( ) );
  
  // Find gap in super.
  int super_id = -1;
  int super_pos = -1;
  int gap_size = -1;
  int gap_dev = -1;
  for (int ii=0; ii<supers.isize( ); ii++) {
    const superb &sup = supers[ii];
    for (int jj=0; jj<sup.Ntigs( )-1; jj++) {
      if ( sup.Tig( jj ) == CONTIG_ID ) {
	super_id = ii;
	super_pos = jj;
	gap_size = sup.Gap( jj );
	gap_dev = sup.Dev( jj );
	break;
      }
    }
  }
  if ( super_id < 0 ) {
    cout << "Fatal error: c" << CONTIG_ID << " not found.\n" << endl;
    return;
  }
  int last_tig = supers[super_id].Ntigs( ) - 1;
  int cgleft_id = supers[super_id].Tig( super_pos );
  int cgright_id = supers[super_id].Tig( super_pos + 1 );

  log << "left cg:   c" << cgleft_id
      << " = s" << super_id << "." << super_pos << "/" << last_tig
      << "\t" << supers[super_id].Len( super_pos ) << " bp\n"
      << "gap:       [" << gap_size << " +/- " << gap_dev << "]\n"
      << "right cg:  c" << cgright_id
      << " = s" << super_id << "." << super_pos + 1 << "/" << last_tig
      << "\t" << supers[super_id].Len( super_pos + 1 ) << " bp\n" << endl;
  
  // Load (with SparseRead) KmerPaths of contigs.
  vecKmerPath cg_paths;
  {
    vec<int> select = MkVec( cgleft_id, cgright_id );
    sort( select.begin( ), select.end( ) );
    cg_paths.SparseRead( contigs_paths_file, select, 0 );
  }
  
  // Find kmer ids to walk.
  longlong kid_left = -1;
  tagged_rpint tag_left;
  {
    bool tag_found = false;
    const KmerPath &left_kp = cg_paths[cgleft_id];
    for (int ii=left_kp.NSegments( )-1; ii>=0; ii--) {
      if ( tag_found ) break;
      const KmerPathInterval &kpi = left_kp.Segment( ii );
      for (longlong kid=kpi.Stop( ); kid>=kpi.Start( ); kid--) {
	if ( tag_found ) break;
	vec<longlong> left_locs;
	Contains( unipathsdb, kid, left_locs );
	for (int jj=0; jj<left_locs.isize( ); jj++) {
	  if ( tag_found ) break;
	  if ( unipathsdb[ left_locs[jj] ].Rc( ) ) continue;
	  kid_left = kid;
	  tag_left = unipathsdb[ left_locs[jj] ];
	  tag_found = true;
	}
      }
    }
    ForceAssert( tag_found );
  }

  longlong kid_right = -1;
  tagged_rpint tag_right;
  {
    bool tag_found = false;
    const KmerPath &right_kp = cg_paths[cgright_id];
    for (int ii=0; ii<right_kp.NSegments( ); ii++) {
      if ( tag_found ) break;
      const KmerPathInterval &kpi = right_kp.Segment( ii );
      for (longlong kid=kpi.Start( ); kid<=kpi.Stop( ); kid++) {
	if ( tag_found ) break;
	vec<longlong> right_locs;
	Contains( unipathsdb, kid, right_locs );
	for (int jj=0; jj<right_locs.isize( ); jj++) {
	  if ( tag_found ) break;
	  if ( unipathsdb[ right_locs[jj] ].Rc( ) ) continue;
	  kid_right = kid;
	  tag_right = unipathsdb[ right_locs[jj] ];
	  tag_found = true;
	}
      }
    }
    ForceAssert( tag_found );
  }
  
  // Walk on graph.
  int start = tag_left.ReadId( );
  int target = tag_right.ReadId( );
  
  int level = 0;
  bool target_found = false;
  vec<int> traversed( 1, start );
  vec<int> level_vtx( 1, start );
  while ( level < DEPTH ) {
    level++;
    if ( level_vtx.size( ) < 1 ) break;

    vec<int> next;
    for (int ii=0; ii<level_vtx.isize( ); ii++) {
      const vec<int> &adder = adjacency.From( level_vtx[ii] );
      copy( adder.begin( ), adder.end( ), back_inserter( next ) );
    }
    sort( next.begin( ), next.end( ) );
    next.erase( unique( next.begin( ), next.end( ) ), next.end( ) );
    
    level_vtx.clear( );
    for (int ii=0; ii<next.isize( ); ii++) {
      int vtx = next[ii];
      bool seen = binary_search( traversed.begin( ), traversed.end( ), vtx );
      if ( ! target_found ) target_found = ( target == vtx );
      int uni_cn = EstimatedCN( copynum[vtx] );
      log << "lev" << level << "." << ii << "/" << next.size( ) - 1 << "\t"
	  << "v=" << vtx << "\t"
	  << unipaths_paths[vtx].TotalLength( ) << " kmers, cn="
	  << uni_cn << ( target_found ? "\tTARGET FOUND" : "" );
      if ( seen ) log << "\tLOOP\n";
      if ( uni_cn > MAX_CN ) log << "\tCUT - HIGH CN BRANCH\n";
      if ( ! seen ) {
	traversed.push_back( vtx );
	if ( uni_cn <= MAX_CN ) level_vtx.push_back( vtx );
	log << "\n";
      }
    }
    log << "\n";
    
    if ( target_found ) break;
    if ( level < DEPTH) sort( traversed.begin( ), traversed.end( ) );
  }
  
  // Generate dot, eps.
  if ( EPS ) {
    vec<String> vlabels( adjacency.N( ), "" );
    for (int ii=0; ii<traversed.isize( ); ii++) {
      int vtx = traversed[ii];
      int uni_cn = EstimatedCN( copynum[vtx] );
      String str_vtx = ToString( vtx );
      String str_len = ToString( unipaths_paths[vtx].TotalLength() ) + " kmers";
      String str_cn = "cn: " + ToString( uni_cn );
      vlabels[vtx] =  str_vtx + "\\n" + str_len + "\\n" + str_cn;
      if ( vtx == start ) vlabels[vtx] += "\\nSTART";
      if ( vtx == target ) vlabels[vtx] += "\\nSTOP";
    }    
    ofstream dotout( dot_file.c_str( ) );
    adjacency.DOT_vl( dotout, vlabels, traversed );
    dotout.close( );

    String theCommand = "dot -Tps " + dot_file + " -o " + eps_file;
    RunCommand( theCommand );
    
    theCommand = "gv " + eps_file + " &";
    RunCommand( theCommand );
  }

  // Done.
  log.close( );
  
}

