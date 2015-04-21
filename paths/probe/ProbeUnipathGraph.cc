///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// paths/probe/ProbeUnipathGraph.cc
//
// Patch the gaps in unipath scaffolds using the global unipath graph.
// This module is similar to KPatch, but we require all the contigs in
// the scaffolds to be ( unique? ) unibases.
//
// We align one end of the jump read to the scaffold and the other side
// points to unipaths that can potentially be used to patch the gap. Using 
// these constraints to find the suitable paths that have the correct
// size to fill the gap.
//
// To simplify the data structure, all the jump reads are aligned to the
// unibases in only one direction. If we consider a pair of jump reads "starts"
// from a position x1 at unibase u1, and "stops" at at position x2 at unibase
// u2, what we really care is the value of len1 - x1 and x2. We can record 
// the alignments as r2 aligned to u2 at position x2, and r1 aligned to u1*
// at position len1 - x1, both in the the same direction.
//
// As show in the diagram below, we are converting
//                         -->@              @<-- 
// ----------u1--------------------       ------------------u2-------------
//
// to
//     @<--                                  @<-- 
// ----------u1*-------------------       ------------------u2-------------
//
// In addition, we are only intersted in the position of the "high quality" end
// of the of reads ( marked as @ in the diagram for shear jumps), which should
// be better defined compared with the low quality end alignment position.
//
// * Not the the arrow of the reads represents the direction that the reads
// are written in the library using allpaths convention as "innies".
//

// MakeDepend: library OMP
#include <omp.h>
#include <unistd.h>
#include "MainTools.h"
#include "Basevector.h"
#include "FastIfstream.h"
#include "PairsManager.h"
#include "Superb.h"
#include "kmers/KmerRecord.h"
#include "math/Functions.h"
#include "lookup/FirstLookupFinderECJ.h"
#include "paths/PdfEntry.h"
#include "paths/GetNexts.h"
#include "paths/UnibaseUtils.h"
#include "paths/reporting/ReftigUtils.h"
#include "ParallelVecUtilities.h"
#include "paths/SaveScaffoldGraph.h"


int main(int argc, char *argv[])
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_String_OrDefault(HEAD, "extended40.shaved.patched");
  CommandArgument_Int_OrDefault(K, 96);
  CommandArgument_Int_OrDefault(LOOKUP_K, 12);
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
      "Number of threads to use (use all available processors if set to 0)");
  CommandArgument_String_OrDefault_Doc(TIGS, "",
      "if specified, only process gaps have these contigs on their left");
  CommandArgument_String_OrDefault(JUMP_READS, "jump_reads_filt_cpd");
  CommandArgument_Bool_OrDefault(WRITE, True);
  // FLIP = True for shear-jump libraries
  CommandArgument_Bool_OrDefault_Doc(FLIP, True,
    "If true, rc reads, error correct, and then rc back" );

  // Logging

  CommandArgument_Int_OrDefault_Doc(VERBOSITY, -1,
      "VERBOSITY=1: print one line for each gap\n"
      "VERBOSITY=2: print more data for each gap\n"
      "VERBOSITY=3: print a lot more\n"
      "VERBOSITY=4: print even more\n"
      "(Note VERBOSITY>=2 should only be used if TIGS is set to one contig.)");
  CommandArgument_String_OrDefault(DEBUG_READS, "");
  CommandArgument_Bool_OrDefault(VALIDATE, False);
  CommandArgument_Bool_OrDefault_Doc(PERFECT, False,
      "Only display perfect alignment to reference");
  CommandArgument_Bool_OrDefault_Doc(VIS_ALIGN, False,
      "Display unibases alignment to reference");
  CommandArgument_String_OrDefault_Doc(GRAPH_OUT, "",
      "Output the unipath graphs (mapped to reference ) ");

  EndCommandArguments;

  // Define directories

  cout << Date( ) << ": begin" << endl;
  String data_dir = PRE + "/" + DATA;
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String tmp_dir = run_dir + "/tmp";

  // input files

  // jump_reads
  String jquals_fn = run_dir + "/" + JUMP_READS + ".qualb";
  String jreads_fn = run_dir + "/" + JUMP_READS + ".fastb";
  String jpairs_fn = run_dir + "/" + JUMP_READS + ".pairs";
  // unibases 
  String KS = ToString(K);
  String unibase_fn    = run_dir + "/" + HEAD + ".unibases.k" + KS ;
  String unibase_CN_fn = run_dir + "/" + HEAD + ".unipaths.predicted_count.k" + KS;
  String unibase_lk_fn = unibase_fn + ".lookup";

  // Thread control

  NUM_THREADS = configNumThreads(NUM_THREADS);
  omp_set_num_threads( NUM_THREADS );

  // Define unibases to process.

  vec<int> tigs_to_process;
  if ( TIGS != "" ) {   
    ParseIntSet( TIGS, tigs_to_process );
    if ( VERBOSITY < 0 ) VERBOSITY = 1;    
  }
  if ( VERBOSITY >= 2 && !tigs_to_process.solo( ) ) {
    cout << "If you use VERBOSITY >= 2, then you need to specify a single contig.\n";
    exit(1);    
  }

  // Load unibases and ancillary data structures.

  cout << Date() << ": loading unibases from " << unibase_fn << endl;
  vecbasevector unibases( unibase_fn );
  int nuni = unibases.size();
  cout << Date() << ":     " << nuni << " unibases loaded." << endl;
  vec<int> to_rc;
  UnibaseInvolution( unibases, to_rc );
  vec< vec<int> > nexts;
  GetNexts( K, unibases, nexts );

  // Get copy number.

  //cout << Date( ) << ": getting copy number" << endl;
  //VecPdfEntryVec CN( unibase_CN_fn.c_str( ) );
  //vec<int> predicted_CN( nuni, -1 );
  //for ( int i = 0; i < nuni; i++ ) GetMostLikelyValue( predicted_CN[i], CN[i] );
  //vec<double> CN_raw;
  //BinaryReader::readFile( run_dir + "/" + HEAD + ".unipaths.cn_raw.k" + KS, &CN_raw );
     
  // Align the jumping reads to these unibases.  
  // Assuming allpaths convention is used ( reads are stored as innies ), 
  // it's better to flip the shear jump reads so that the high quality end is at the 
  // beginning of the reads and record the alignment position from this end

  PairsManager pairs( jpairs_fn );
  cout << Date() << ": loading jump reads from " << jreads_fn << endl;
  vecbasevector jreads( jreads_fn );
  size_t njreads = jreads.size();
  cout << Date() << ":     " << njreads << " jump reads loaded." << endl;

  cout << Date() << ": loading jump reads quals from " << jquals_fn << endl;
  VecQualNibbleVec jquals;
  LoadQualNibbleVec(jquals_fn, &jquals);
  ForceAssertEq(njreads, jquals.size());

  // Flip the bases and quals if necessary.
  if (FLIP) {
    cout << Date() << ": rc-ing query bases" << endl;
    for (size_t ii = 0; ii < njreads; ii++) jreads[ii].ReverseComplement();
    cout << Date() << ": rc-ing query quals" << endl;
    for (size_t ii = 0; ii < njreads; ii++) jquals[ii].ReverseMe();
  }    
  
  // Run FirstLookup.
  vec<first_look_align> first_aligns;
  FirstLookupFilterECJ lookup_filter;
  {
    // align the reads using FirstLookAlignECJ
    int MIN_MATCH=20, MISMATCH_THRESHOLD=3, MISMATCH_NEIGHBORHOOD=8, MISMATCH_BACKOFF=3, SCORE_DELTA=20, SCORE_MAX=100, MAX_KMER_FREQ=10000;
    // Set up filtering options for call to FirstLookupFinedECJ.
    //lookup_filter.orientation = FLIP ?  FirstLookupFilterECJ::RC_ONLY : FirstLookupFilterECJ::FW_ONLY;
    lookup_filter.orientation = FirstLookupFilterECJ::FW_ONLY; // always forward aligned to the unibases
    lookup_filter.min_size = 20;
    lookup_filter.max_kmer_freq = MAX_KMER_FREQ;
    lookup_filter.max_extend = MAX_KMER_FREQ;
    lookup_filter.score_delta = SCORE_DELTA;
    lookup_filter.score_max = SCORE_MAX;
    lookup_filter.min_match = MIN_MATCH;
    lookup_filter.mismatch_threshhold = MISMATCH_THRESHOLD;
    lookup_filter.mismatch_neighborhood = MISMATCH_NEIGHBORHOOD;
    lookup_filter.mismatch_backoff = MISMATCH_BACKOFF;
    lookup_filter.max_placements = 2;   // if there are multiple placements
    ParseIntSet(DEBUG_READS, lookup_filter.debug_reads, False);
  }

  {
    // Load LookupTab for alignment.
    cout << Date() << ": Loading LookupTab from file " << unibase_lk_fn << endl;
    LookupTab lookup_tab(unibase_lk_fn.c_str());
    ForceAssertEq( (unsigned int )LOOKUP_K, lookup_tab.getK() );
    cout << Date() << ": running FirstLookupFinderECJ" << endl;
    FirstLookupFinderECJ lfinder(lookup_filter, lookup_tab, unibases, K); // I don't think K here will make any difference
    lfinder.getAllAlignments(jreads, jquals, &first_aligns, NUM_THREADS);
  }
  cout << Date() << ":     " << first_aligns.size() << " alignments found from jump reads to unibases." << endl;

  // Load reference genome for validation purpose

  if ( VALIDATE ) {
    String genome_fn        = data_dir + "/genome.fastb";
    String genome_lookup_fn = data_dir + "/genome.lookup"; // the LookupTable
    vecbasevector genome( genome_fn );
    // align reads to reference
    //String genome_lktab_fn  = run_dir + "/genome.lookup";  // the LookupTab  for FirstLookupFinderECJ
    //vec<first_look_align> jump2genome;
    //{
    //  // Load LookupTab for alignment.
    //  cout << Date() << ": Loading LookupTab from file " << genome_lktab_fn<< endl;
    //  LookupTab lookup_tab(genome_lktab_fn.c_str());
    //  ForceAssertEq( (unsigned int )LOOKUP_K, lookup_tab.getK() );
    //  cout << Date() << ": running FirstLookupFinderECJ" << endl;
    //  FirstLookupFinderECJ lfinder(lookup_filter, lookup_tab, genome, K);
    //  lfinder.getAllAlignments(jreads, jquals, &jump2genome, NUM_THREADS);
    //}
    //cout << Date() << ":     " << jump2genome.size() << " alignments found from jump reads to reference." << endl;

    // align unibases to reference genome

    String aligns_file = tmp_dir + "/unibases_to_genome.qltout";
    vec<look_align> unibase_aligns;
    GetAlignsFast( K, unibase_fn, genome_lookup_fn, aligns_file, unibase_aligns, True, tmp_dir);
    cout << Date() << ":     " << unibase_aligns.size() << " alignments found from unibases to reference." << endl;

    // Remove RC alignment to genome
    vec<Bool> toErase( unibase_aligns.size(), False );
    for( int i = 0; i < unibase_aligns.isize(); i++ ) 
      if ( ( PERFECT && unibase_aligns[i].Errors() != 0 ) || unibase_aligns[i].IsQueryRC() ) 
	toErase[i] = True;
    EraseIf( unibase_aligns, toErase );
    order_lookalign_Begin sorter;
    ParallelSort( unibase_aligns, sorter );

    // making a fake scaffold

    vec< vec<String> > table; // for scaffold visualization
    cout << Date() << ": Generating fake scaffolds. " << endl;
    vec< superb > fake_scaffolds;
    int prev_target_id = -1;
    for (size_t ii=0; ii< unibase_aligns.size( ); ii++) {
      if ( unibase_aligns[ii].IsQueryRC( ) ) continue;
      int target_id = unibase_aligns[ii].target_id;
      int query_id =  unibase_aligns[ii].query_id;

      // The gap table
      if ( target_id != prev_target_id ) { // new scaffold
	fake_scaffolds.push_back( superb() );
	fake_scaffolds.back().PlaceFirstTig( query_id, unibases[query_id].size() );
	if ( ! table.empty() ) {
	  cout << "----scaffold aligned to genome id " << prev_target_id << endl;
	  if ( VERBOSITY >= 2 )
	       PrintTabular( cout, table, 1, "lrlll" );
	}
	table.clear();
	table.push_back( MkVec( ToString(query_id), ToString(unibases[query_id].size()), 
		       String("-"), ToString(unibase_aligns[ii].pos2()), ToString(unibase_aligns[ii].Pos2()) ) );
	prev_target_id = target_id;
	continue;
      } 
      superb& the_scaffold = fake_scaffolds.back();
      int gap = unibase_aligns[ii].pos2() -  unibase_aligns[ii-1].Pos2(); 
      the_scaffold.AppendTig(query_id, unibases[query_id].size(),  gap, 0);
      table.push_back( MkVec( ToString(query_id), ToString(unibases[query_id].size()), 
		     ToString(gap), ToString(unibase_aligns[ii].pos2()), ToString(unibase_aligns[ii].Pos2()) ) );
    }
    cout << Date( ) << ": saving supers" << endl;
    String supers_file = "fake_scaffolds.superb";
    WriteSuperbs( supers_file, fake_scaffolds);
    cout << Date() << ":     " << fake_scaffolds.size() << " scaffolds generated from reference." << endl;

    // Visualization of Unibases on reference
    if ( VIS_ALIGN ) { 
      cout << Date() << ": visualize unibases on reference. " << endl;
      vec< vec<String> > vis_aligns_all( genome.size() );  // one table for each genome id
      for( size_t i = 0; i < vis_aligns_all.size(); i++ )
	vis_aligns_all[i].push_back( genome[i].ToString() );
      // The table
      for (size_t ii=0; ii< unibase_aligns.size( ); ii++) {
	look_align& a = unibase_aligns[ii];
	if ( a.IsQueryRC( ) ) continue;
	int target_id = a.target_id;
	int query_id =  a.query_id;
	size_t start =  a.pos2();
	size_t end   =  a.Pos2();

	const basevector& query = unibases[ query_id ];
	const basevector& target = genome[ target_id ];
	int shift = a.pos1();
	String l1; // alignment sequence to be displayed
	for ( int k = 0; k < shift; k++ ) l1.push_back( as_base( query[k] ) );
	String l1_to_add; // additional info to be displayed
	int p1 = a.pos1();
	int p2 = a.pos2();
	for ( int j = 0; j < a.a.Nblocks(); ++j ){
	  int gap = a.a.Gaps(j);
	  int len = a.a.Lengths(j);
	  if ( gap > 0 ) {
	    l1.append( gap, 'D' );
	  }
	  if ( gap < 0 ) {
	    int nchar = -gap;
	    l1_to_add.push_back( '^' );
	    for ( int x = 0; x < nchar; x++ )
	      l1_to_add.push_back( as_base( query[p1++]) );
	    l1_to_add.push_back( '$' );
	  }
	  // are the perfect matching blocks long enough for outputing previous insertion 
	  int avail_counter = 0;
	  int avail_start = 0;
	  if ( len > (int) l1_to_add.size() ) {
	    for ( int x = 0, counter = 0, p= 0; x < len; x++ ) 
	      if ( query[p1+x] != target[p2+x] 
		  || x == 0 && gap < 0 ) { // reset, also exclude first char after reference gap
		counter = 0;
		p= x+1;
	      } else{
		counter++;
		if ( counter > avail_counter ) { 
		  avail_counter = counter; 
		  avail_start = p;
		}
	      }
	  }
	  // output the matching block
	  for ( int x = 0; x < len; x++, p1++, p2++ ) {
	    //l2.push_back( as_base( target[p2] ) );
	    if ( x == 0 && gap < 0 ) {
	      l1.push_back( tolower( as_base( query[p1] ) ) );
	    }else if ( query[p1] != target[p2] ) {
	      l1.push_back( as_base( query[p1] ) );
	    }else { // in the matching block
	      if ( l1_to_add.size() > 0 && avail_counter >= (int) l1_to_add.size() 
		  && x >= avail_start && x < avail_start + (int) l1_to_add.size() ) {
		l1.push_back( l1_to_add[x - avail_start] );
	      }else {
		l1.push_back( ( '-' ) );
	      }
	    }
	  }
	  if ( avail_counter >= (int) l1_to_add.size() ) l1_to_add.clear();
	}
	// the connectivity
	l1.push_back( '#' );
	for ( int j = 0; j < nexts[query_id].isize(); j++ ) 
	  l1.append( ToString( nexts[query_id][j] ) + "," );
	if ( l1_to_add.size() > 0 )  // add the rest of inserts to the end
	  l1.append( l1_to_add );

	String header;
	// if the first K-1 mer aligns to ref perfectly, we can ignore them to simplify the diagram
	bool toIgnore = true;
	for ( int j = 0; j < K - 1; j++ ) if ( l1[j] != '-' ) { toIgnore = false; break; }
	if ( toIgnore ) {
	  header  = "<<" + ToString(query_id) + "#";
	  l1 = l1.substr( K-1, l1.size() - K + 1 );
	  start += K - 1;
	}else{
	  header = ToString(query_id) + "#";
	}
	if ( query.size() > 500 ) header.append( "L" + ToString( query.size() ) );
	shift += header.size();

	// find a line where the total length is smaller than pos2 and append
	// add new line if necessory
	vec<String> & vis_aligns = vis_aligns_all[target_id];
	size_t line_index = 0; 
	for ( ; line_index < vis_aligns.size(); ++line_index )
	  if ( vis_aligns[line_index].size() + shift < start ) break;
	if ( line_index == vis_aligns.size() ) vis_aligns.push_back(String());

	ForceAssertGe( start, shift + vis_aligns[line_index].size() );
	vis_aligns[line_index].append( String( start - shift - vis_aligns[line_index].size(), ' ' ) );

	vis_aligns[line_index].append(header);
	vis_aligns[line_index].append(l1);
      }
      // print the table
      for ( size_t i = 0; i < vis_aligns_all.size(); ++i ) {
	cout << "-------------- alignment to genome " << i << endl;
	for ( size_t j = 0; j < vis_aligns_all[i].size(); ++j )
	  cout << vis_aligns_all[i][j] << endl;
      }
    }
    // check the unipath connection
    vec< triple<int64_t,int,int> > end_pos( unibases.size(), triple<int64_t,int,int>(-1,-1,-1) );
    for ( size_t ii=0; ii< unibase_aligns.size( ); ii++ ) {
      look_align& a = unibase_aligns[ii];
      int target_id = a.target_id;
      int query_id =  a.query_id;
      for ( int j = 0; j < nexts[query_id].isize(); j++ ) 
	end_pos[ nexts[query_id][j] ] = triple<int64_t,int,int>(a.Pos2() - K + 1, target_id, query_id);
      if ( end_pos[query_id].first > 0 &&  end_pos[query_id].first != a.pos2() ) {
	cout << "Connection missmatch from " << end_pos[query_id].third << " to " << query_id << endl;
	cout << "   starting pos " << a.pos2() << " at target " << target_id << endl;
        cout << "   expected pos " <<  end_pos[query_id].first << " at target " << end_pos[query_id].second <<  endl;
      }
    }

    // output the graph ( of the aligned unibases only )

    if ( GRAPH_OUT != "" )
    {
      ofstream out( GRAPH_OUT.c_str() );
      // only take the edges that are aligned to reference
      vec <int> newindexs;
      for ( size_t i = 0; i < unibase_aligns.size(); ++i ) {
	newindexs.push_back( unibase_aligns[i].query_id );
      }
      UniqueSort( newindexs );
      vec <int> mapping( nuni, -1 );
      for ( size_t i = 0; i < newindexs.size(); ++i ) 
	mapping[ newindexs[i] ] = i;
      cout << "Reference-aligned unibases = " << newindexs.size() << endl;
      // from and to in the remaining unibases
      size_t nuni2 = newindexs.size();
      vec< vec<int> > from ( nuni2 ); 
      vec< vec<int> > to ( nuni2 );
      for ( size_t i = 0; i < nuni2; ++i ) {
	int iold = newindexs[i];
	for ( size_t j = 0; j < nexts[ iold ].size(); ++j ) {
	  int jnew = mapping[ nexts[iold][j] ];
	  if ( jnew < 0 ) continue;
	  to[ jnew ].push_back( i );
	  from[ i ].push_back( jnew );
	}
      }
      // prepare the graph that use unibases as edges
      vec <basevector>  edges( nuni2 );
      for ( size_t i = 0; i < unibases.size(); ++i ) 
	if ( mapping[i] >= 0 ) edges[ mapping[i] ] = unibases[i];
      // edge_index[2*i] and [2*i+1] is the start/end of the edge
      vec <int> edge_index( edges.size() * 2, vec <int> :: IDENTITY ); 
      vec< vec<int> > groups( edges.size() * 2 );
      for ( size_t i = 0; i < groups.size(); ++i ) groups[i].push_back(i);
      for ( size_t i = 0; i < to.size(); ++i ) {
	for ( size_t j = 0; j < to[i].size(); ++j ) {
	  int v = to[i][j]; // merge all v_end and i_begin nodes
	  int group_id_i = edge_index[ 2 * i ];
	  int group_id_v = edge_index[ 2 * v + 1 ];
	  for ( size_t k = 0; k < groups[ group_id_v ].size(); ++k ) {
	    edge_index[  groups[ group_id_v ][k] ] = group_id_i;
	  }
	  groups[ group_id_i ].append( groups[group_id_v] );
	  groups[ group_id_v ].clear();
	}
      }
      // delete empty groups, and adjust index
      vec <Bool> toDel( groups.size(), False );
      for ( size_t i = 0; i < groups.size(); ++i ) toDel[i] = groups[i].empty();
      EraseIf( groups, toDel );
      for ( size_t i = 0; i < groups.size(); ++i ) 
	for ( size_t j = 0; j < groups[i].size(); ++j ) 
	  edge_index[ groups[i][j] ] = i;
      // the vertices
      vec <int> joints;
      for ( size_t i = 0; i < edge_index.size(); ++i ) 
	if ( edge_index[i] >= 0 ) joints.push_back( edge_index[i] );
      UniqueSort( joints );
      unsigned nvertex = joints.size();
      cout << "Number of vertices is " << nvertex << endl;
      // prepare the graph for outputing
      vec< vec<int> > new_to( nvertex );
      vec< vec<int> > new_from( nvertex );
      vec< vec<int> > to_edge_obj( nvertex );
      vec< vec<int> > from_edge_obj( nvertex );
      // unibases as edges
      for ( size_t i = 0; i < to.size(); ++i ) {
	int start_index = edge_index[ 2 * i ];
	int end_index = edge_index[ 2 * i + 1 ];
	ForceAssertGe( start_index, 0 );
	ForceAssertGe( end_index, 0 );
	new_from[ start_index ].push_back( end_index );
	from_edge_obj[ start_index ].push_back( i );
	new_to[ end_index ].push_back( start_index );
	to_edge_obj[ end_index ].push_back( i );
      }

      for ( size_t i = 0; i < new_to.size(); ++i ) {
	SortSync( new_to[i], to_edge_obj[i] );
	SortSync( new_from[i], from_edge_obj[i] );
      }
      digraphE <basevector> graph( new_from, new_to, edges, to_edge_obj, from_edge_obj );
      Bool label_contigs = True;
      Bool label_vertices = True;
      vec <double> lengths( edges.size(), 0);
      for ( size_t i = 0; i < edges.size(); ++i ) lengths[i] = edges[i].size();
      graph.PrettyDOT( out, lengths, digraphE<basevector>::edge_label_info( 
          digraphE<basevector>::edge_label_info::EXPLICIT, False, False, NULL),
          label_contigs, label_vertices );


      out.close();
    }
  }

  cout << "\ncommand: " << command.TheCommand( ) << "\n\n";    
}
