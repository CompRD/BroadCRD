///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Bitvector.h"
#include "CoverageAnalyzer.h"
#include "SeqInterval.h"
#include "graph/Digraph.h"
#include "math/NStatsTools.h"
#include "paths/Unipath.h"
#include "util/SearchFastb2Core.h"
#include "util/StroboAligner.h"
#include "feudal/BinaryStream.h"

void PrintStatsCore( ostream &out,
		     const int K,
		     const int super_id,
		     const int super_klen,
		     const String &descriptor_base,
		     const CoverageAnalyzer &analyzer,
		     const digraphE<int> *graphE = 0,
		     const vec<int> *to_left = 0,
		     const vec<int> *to_right = 0,
		     const vec<seq_interval> *super_gaps = 0,
		     const vec< pair<int,int64_t> > *start2id = 0 );

/**
 * EvalUnipahtsPerfect
 *
 * Eval unipaths by aligning them onto a reference. This is a two
 * steps procedure: first unibases are aligned onto a reference
 * supercontig, and then gaps on the reference are filled by walking
 * on the graph.
 *
 * K: unipaths kmer size
 * REF_HEAD: it will load <REF_HEAD>.{fastb,fastamb}
 * REF_ID: if of selected super
 * QUERY_HEAD: it needs various <QUERY_HEAD>.{*}.k<K> files
 * OUT_DIR: full path name of output dir
 * ML: ie MAX_LOOPS, argument of AllPathsLengthRange in digraphE
 * MIN_LENGTH: in kmers, min length of unibases from which gaps are bridged
 * STROBE*: StroboAligner arguments
 * OVERWRITE: overwrite existing files
 * VERBOSE: log's verbose flag
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_Int( K );
  CommandArgument_String( REF_HEAD );
  CommandArgument_Int( REF_ID );
  CommandArgument_String( QUERY_HEAD );
  CommandArgument_String( OUT_DIR );
  CommandArgument_Int_OrDefault( ML, 10000 );
  CommandArgument_Int_OrDefault( MIN_LENGTH, 64 );
  CommandArgument_Int_OrDefault( STROBE_K, 96 );
  CommandArgument_Int_OrDefault( STROBE_INTERVAL, 100 );
  CommandArgument_Int_OrDefault( STROBE_MAX_GAP, 40 );
  CommandArgument_Int_OrDefault( STROBE_MIN_LENGTH, 1000 );
  CommandArgument_Double_OrDefault( STROBE_MIN_RATIO, 0.8 );
  CommandArgument_Bool_OrDefault( OVERWRITE, False );
  CommandArgument_Bool_OrDefault( VERBOSE, False );
  EndCommandArguments;
  
  // File names.
  String strK = "k" + ToString( K );
  
  String target_bases_file = REF_HEAD + ".fastb";
  String target_amb_file = REF_HEAD + ".fastamb";
  String query_hkp_file = QUERY_HEAD + ".hyper." + strK;
  String query_hkpint_file = QUERY_HEAD + ".hyper_int." + strK;
  String query_bases_file = QUERY_HEAD + ".unibases." + strK;
  String query_unipaths_file = QUERY_HEAD + ".unipaths." + strK;
  String query_unipathsdb_file = QUERY_HEAD + ".unipathsdb." + strK;
  String query_paths_file = QUERY_HEAD + ".paths." + strK;
  String query_paths_rc_file = QUERY_HEAD + ".paths_rc." + strK;
  String query_pathsdb_file = QUERY_HEAD + ".pathsdb." + strK;
  String query_adj_file = QUERY_HEAD + ".unipath_adjgraph." + strK;

  String log_file = OUT_DIR + "/EvalUnipathsPerfect.log";
  String select_target_bases_file = OUT_DIR + "/target.fastb";
  String strobo_aligns_file = OUT_DIR + "/strobo_aligns.out";
  String aligns_file = OUT_DIR + "/aligns.out";

  Mkpath( OUT_DIR );

  ofstream log( log_file.c_str( ) ); 
  PrintCommandPretty( log );
  cout << "Sending log to " << log_file << endl;

  // Load reference (NOTICE: super lengths and gaps in kmers!)
  vec<int> super_lens;
  vec<seq_interval> super_gaps;
  {
    log << Date( ) << ": loading super gaps" << endl;
    vecbitvector amb( target_amb_file.c_str( ) );
    
    super_lens.reserve( amb.size( ) );
    for (int ii=0; ii<(int)amb.size( ); ii++)
      super_lens.push_back( (int)amb[ii].size( ) - ( K - 1 ) );

    for (int super_id=0; super_id<(int)amb.size( ); super_id++) {
      int pos = 0;
      int cgpos = 0;
      while ( pos < (int)amb[super_id].size( ) ) {
	while( pos < (int)amb[super_id].size( ) && ! amb[super_id][pos] ) pos++;
	int begin = Max( 0, pos - ( K - 1 ) );
	while( pos < (int)amb[super_id].size( ) && amb[super_id][pos] ) pos++;
	int end = pos;
	if ( begin < (int)amb[super_id].size( ) ) {
	  seq_interval sigap( cgpos, super_id, begin, end );
	  super_gaps.push_back( sigap );
	  cgpos++;
	}
      }
    }
  }

  // REF_ID length in kmers.
  const int slen = super_lens[REF_ID];
  
  // Load query.
  digraphE<int> graphE;
  
  if ( IsRegularFile( query_hkpint_file ) ) {
    log << Date( ) << ": loading digraphE<int>" << endl;
    BinaryReader::readFile( query_hkpint_file, &graphE );
  }
  else {
    HyperKmerPath hkp;
    
    if ( IsRegularFile( query_hkp_file ) ) {
      log << Date( ) << ": loading hyper kmer path" << endl;
      BinaryReader::readFile( query_hkp_file, &hkp );
    }
    else {
      log << Date( ) << ": loading unipaths" << endl;
      vecKmerPath unipaths( query_unipaths_file );
      
      digraph graph; 
      if ( ! IsRegularFile( query_adj_file ) ) {
	log << Date( ) << ": loading paths data to build adj graph" << endl;
	vecKmerPath paths( query_paths_file );
	vecKmerPath paths_rc( query_paths_rc_file );    
	BREAD2( query_pathsdb_file, vec<tagged_rpint>, pathsdb );

	log << Date( ) << ": loading unipathsdb to build adj graph" << endl;
	BREAD2( query_unipathsdb_file, vec<tagged_rpint>, unipathsdb );
	
	log << Date( ) << ": building adjaceny graph" << endl;
	BuildUnipathAdjacencyGraph( paths, paths_rc, pathsdb,
				    unipaths, unipathsdb, graph );
	
	log << Date( ) << ": saving adjaceny graph" << endl;
	BinaryWriter::writeFile( query_adj_file, graph );
      }
      else {
	log << Date( ) << ": loading adjacency graph" << endl;
	BinaryReader::readFile( query_adj_file, &graph );
      }

      log << Date( ) << ": turning digraph into digraphE" << endl;
      BuildUnipathAdjacencyHyperKmerPath( K, graph, unipaths, hkp );
      
      log << Date( ) << ": saving digraphE" << endl;
      BinaryWriter::writeFile( query_hkp_file, hkp );
    }
    
    log << Date( ) << ": generating digraphE<int>" << endl;
    vec<int> k_lens;
    k_lens.reserve( hkp.EdgeObjectCount( ) );
    for(int ii=0; ii<hkp.EdgeObjectCount( ); ii++)
      k_lens.push_back( hkp.EdgeObject( ii ).TotalLength( ) );
    graphE.Initialize( hkp, k_lens );

    log << Date( ) << ": saving digraphE<int>" << endl;
    BinaryWriter::writeFile( query_hkpint_file, graphE );
  }
  
  // Save selected contigs from reference.
  if ( OVERWRITE || ( ! IsRegularFile( select_target_bases_file ) ) ) {
    log << Date( ) << ": loading target sequence" << endl;
    vecbvec t_bases;
    vec<int> select( 1, REF_ID );
    t_bases.SparseRead( target_bases_file, select, 0 );
    
    log << Date( ) << ": saving selected subset of target sequence" << endl;
    t_bases.WriteAll( select_target_bases_file );
  }

  // Generate and/or load strobo-aligns.
  vec< triple<int64_t,int64_t,int> > strobo_aligns;
  if ( OVERWRITE || ( ! IsRegularFile( strobo_aligns_file ) ) ) {
    vec< triple<int64_t,int64_t,int> > saligns;
    StroboAligner( saligns, select_target_bases_file, query_bases_file,
		   OUT_DIR, STROBE_K, STROBE_INTERVAL, STROBE_MAX_GAP,
		   STROBE_MIN_LENGTH, STROBE_MIN_RATIO );
    
    strobo_aligns.reserve( saligns.size( ) );
    
    // Warning: length is expressed in kmer space, not in bases!
    log << Date( ) << ": saving strobo alignments" << endl;
    ofstream alout( strobo_aligns_file.c_str( ) );
    for (size_t ii=0; ii<saligns.size( ); ii++) {
      if ( saligns[ii].third < 0 ) continue;
      strobo_aligns.push_back( saligns[ii] );
      alout << "ALIGN " << saligns[ii].first
	    << " "  << saligns[ii].second
	    << " + " << saligns[ii].third
	    << " " << saligns[ii].third + graphE.EdgeObject(saligns[ii].first)
	    << "\n";
    }
    alout.close( );
  }
  else {
    log << Date( ) << ": loading existing strobo alignments" << endl;
    String keyALIGN, orientation;
    int64_t q_id, t_id;
    int begin, end;
    ifstream alin( strobo_aligns_file.c_str( ) );
    while ( 1 ) {
      alin >> keyALIGN >> q_id >> t_id >> orientation >> begin >> end;
      if ( !alin ) break;
      if ( orientation == "-" ) continue;   // this should not happen
      strobo_aligns.push_back( make_triple( q_id, t_id, begin ) );
    }
    alin.close( );
  }

  // Generate and/or load alignments (will keep only the fw ones).
  vec< triple<int64_t,int64_t,int> > fw_aligns;
  if ( OVERWRITE || ( ! IsRegularFile( aligns_file ) ) ) {
    log << Date( ) << ": starting to align" << endl;
    vec< triple<int64_t,int64_t,int> > aligns;
    SearchFastb2( query_bases_file, select_target_bases_file, K, &aligns );
    fw_aligns.reserve( aligns.size( ) / 2 );
    
    // Warning: length is expressed in kmer space, not in bases!
    log << Date( ) << ": saving alignments (fw only!)" << endl;
    ofstream alout( aligns_file.c_str( ) );
    for (size_t ii=0; ii<aligns.size( ); ii++) {
      if ( aligns[ii].third < 0 ) continue;
      fw_aligns.push_back( aligns[ii] );
      alout << "ALIGN " << aligns[ii].first
	    << " "  << aligns[ii].second
	    << " + " << aligns[ii].third
	    << " " << aligns[ii].third + graphE.EdgeObject( aligns[ii].first )
	    << "\n";
    }
    alout.close( );
  }
  else {
    log << Date( ) << ": loading existing alignments" << endl;
    String keyALIGN, orientation;
    int64_t q_id, t_id;
    int begin, end;
    ifstream alin( aligns_file.c_str( ) );
    while ( 1 ) {
      alin >> keyALIGN >> q_id >> t_id >> orientation >> begin >> end;
      if ( !alin ) break;
      if ( orientation == "-" ) continue;   // this should not happen
      fw_aligns.push_back( make_triple( q_id, t_id, begin ) );
    }
    alin.close( );
  }
  
  // Merge strobo_aligns and fw_aligns.
  log << Date( ) << ": merging aligns, and removing duplicates" << endl;
  vec< triple<int64_t,int64_t,int> > aligns = strobo_aligns;
  copy( fw_aligns.begin( ), fw_aligns.end( ), back_inserter( aligns ) );
  sort( aligns.begin( ), aligns.end( ) );
  aligns.erase( unique( aligns.begin( ), aligns.end( ) ), aligns.end( ) );
  
  // No aligns?
  if ( aligns.size( ) < 1 ) {
    log << "\nNo aligns found.\n" << endl;
    return 0;
  }
  
  // Sort aligns.
  vec< pair<int,int64_t> > start2id;
  start2id.reserve( aligns.size( ) );
  for (size_t align_id=0; align_id<aligns.size( ); align_id++) {
    const triple<int64_t,int64_t,int> &al = aligns[align_id];
    start2id.push_back( make_pair( al.third, al.first ) );
  }
  sort( start2id.begin( ), start2id.end( ) );
  
  // Print statistics of perfect alignments.
  vec<seq_interval> perfect_wins;
  {
    perfect_wins.reserve( start2id.size( ) );
    for (int ii=0; ii<start2id.isize( ); ii++) {
      int query_id = start2id[ii].second;
      int seq_id = REF_ID;
      int begin = start2id[ii].first;
      int end = begin + graphE.EdgeObject( query_id );

      perfect_wins.push_back( seq_interval( query_id, seq_id, begin, end ) );
    }
    
    CoverageAnalyzer analyzer( perfect_wins, &super_lens );
    log << "\n";
    PrintStatsCore( log, K, REF_ID, slen, "PERFECT" , analyzer );
  }
  
  // Merge overlapping unibases.
  vec<int> pillars( 1, 0 );
  for (int ii=start2id.isize( )-1; ii>0; ii--) {
    int current_begin = start2id[ii].first;
    int previous_id = start2id[ii-1].second;
    int previous_begin = start2id[ii-1].first;
    int previous_end = previous_begin + graphE.EdgeObject( previous_id );

    if ( current_begin > previous_end ) pillars.push_back( ii );
    
  }
  sort( pillars.begin( ), pillars.end( ) );

  // All filling paths.
  vec< vec<int> > paths;

  // The new closures.
  vec<seq_interval> imperfect_wins;
  vec< pair<int,int64_t> > start2id_cl;

  // Maps.
  vec<int> toL;
  vec<int> toR;
  graphE.ToLeft( toL );
  graphE.ToRight( toR );

  // Try to fill gaps.
  log << endl;
  for (int pid=0; pid<pillars.isize( )-1; pid++) {

    // Select unibases used as pillars to bridge gap (long ones).
    int s2i_left = -1;
    for (int ii=pillars[pid+1] - 1; ii>=pillars[pid]; ii--) {
      if ( graphE.EdgeObject( start2id[ii].second ) >= MIN_LENGTH ) {
	s2i_left = ii;
	break;
      }
    }
    if ( s2i_left < 0 ) s2i_left = pillars[pid+1] - 1;

    int s2i_right = -1;
    int last = pid<pillars.isize( )-2 ? pillars[pid+2]-1 : start2id.isize( )-1;
    for (int ii=pillars[pid+1]; ii<=last; ii++) {
      if ( graphE.EdgeObject( start2id[ii].second ) >= MIN_LENGTH ) {
	s2i_right = ii;
	break;
      }
    }
    if ( s2i_right < 0 ) s2i_right = pillars[pid+1];

    const pair<int,int64_t> &left = start2id[s2i_left];
    const pair<int,int64_t> &right = start2id[s2i_right];
    
    int64_t left_id = left.second;
    int64_t right_id = right.second;
    int right_begin = right.first;
    int left_begin = left.first;
    int left_end = left_begin + graphE.EdgeObject( left_id );
    int gap_size = right_begin - left_end;
    int v = toR[left_id];
    int w = toL[right_id];
    
    // Info.
    if ( VERBOSE )
      log << "Gap of " << gap_size << " kmers between "
	  << "{v" << toL[left_id] << " e" << left_id << " v" << v << "} "
	  << "len = " << graphE.EdgeObject( left_id ) << " "
	  << "and "
	  << "{v" << w << " e" << right_id << " v" << toR[right_id] << "} "
	  << "len = " << graphE.EdgeObject( right_id ) << "   ";
    
    // Massage min and max allowed gap sizes.
    int minl = gap_size / 5; 
    if ( minl < 10 ) minl = 0;
    int maxl = gap_size + Max( 10, ( gap_size / 5 ) );
    
    // Search for paths (try with Alt method first).
    bool ok = false;
    graphE.AllPathsLengthRangeAlt( v, w, minl, maxl, toR, paths, 1, ML );
    ok = ( paths.size( ) > 0 && paths[0].size( ) > 0 );
    
    if ( !ok ) {
      graphE.AllPathsLengthRange( v, w, minl, maxl, toR, paths, 1, ML );
      ok = ( paths.size( ) > 0 && paths[0].size( ) > 0 );
    }
    
    if ( !ok ) {
      if ( VERBOSE ) log << "NO CLOSURE FOUND\n";
      continue;
    }
    
    // Print info.
    if ( VERBOSE ) {
      log << "--->   [";
      int closure_len = 0;
      for (int ii=0; ii<paths[0].isize( ); ii++) {
	log << paths[0][ii] << ( ii < paths[0].isize( )-1 ? ", " : "]   " );
	closure_len += graphE.EdgeObject( paths[0][ii] );
      }
      log << "closure length = " << closure_len << " kmers\n";
    }
    
    // Add new closure interval.
    seq_interval new_si( paths[0][0], REF_ID, left_end, right_begin );
    imperfect_wins.push_back( new_si );

    int local_len = 0;
    for (int ii=0; ii<paths[0].isize( ); ii++) {
      int64_t edge_id = paths[0][ii];
      int start = left_end + local_len;
      start2id_cl.push_back( make_pair( start, edge_id ) );
      local_len += graphE.EdgeObject( edge_id );
    }
    
  }

  // Add closures to start2ids.
  copy( start2id_cl.begin( ), start2id_cl.end( ), back_inserter( start2id ) );
  sort( start2id.begin( ), start2id.end( ) );
  
  // Stats of perfect and imperfect windows (together).
  vec<seq_interval> wins = perfect_wins;
  copy( imperfect_wins.begin( ), imperfect_wins.end( ), back_inserter(wins) );
  sort( wins.begin( ), wins.end( ) );
  
  CoverageAnalyzer analyzer( wins, &super_lens );
  PrintStatsCore( log, K, REF_ID, slen, "ALL", analyzer,
		  &graphE, &toL, &toR, &super_gaps, &start2id );
  
  // Done.
  log << "\n" << Date( ) << ": done" << endl;
  cout << "\n" << Date( ) << ": done" << endl;
  
}

/**
 * PrintStatsCore
 */
void PrintStatsCore( ostream &out,
		     const int K,
		     const int super_id,
		     const int super_klen,
		     const String &descriptor_base,
		     const CoverageAnalyzer &analyzer,
		     const digraphE<int> *graphE,
		     const vec<int> *to_left,
		     const vec<int> *to_right,
		     const vec<seq_interval> *super_gaps,
		     const vec< pair<int,int64_t> > *start2id )
{
  out << "\nOVERALL STATISTICS FOR " << descriptor_base << " ALIGNMENTS\n";

  vec<seq_interval> cov;
  analyzer.GetCoveragesAtLeast( 1, cov );
  vec<int> cov_lens;
  cov_lens.reserve( cov.size( ) );
  for (int ii=0; ii<cov.isize( ); ii++) {
    if ( cov[ii].SeqId( ) == super_id )
      cov_lens.push_back( cov[ii].Length( ) );
  }
  out << "\n";
  PrintBasicNStats( "covered", cov_lens, out );
  
  vec<seq_interval> gap;
  analyzer.GetCoveragesExactly( 0, gap );
  vec<int> gap_lens;
  gap_lens.reserve( gap.size( ) );
  for (int ii=0; ii<gap.isize( ); ii++) {
    if ( gap[ii].SeqId( ) == super_id )
      gap_lens.push_back( gap[ii].Length( ) );
  }
  out << "\n";
  PrintBasicNStats( "gaps", gap_lens, out );
  
  // Unitigs info.
  if ( ! ( start2id && super_gaps ) ) return;
  
  // All windows (unitigs plus gaps).
  vec<seq_interval> allw;
  allw.reserve( cov.size( ) + gap.size( ) );

  for (int ii=0; ii<cov.isize( ); ii++) {
    seq_interval newsi = cov[ii];
    newsi.SetIntervalId( 1 );
    allw.push_back( newsi );
  }
  for (int ii=0; ii<gap.isize( ); ii++) {
    if ( gap[ii].SeqId( ) != super_id ) continue;
    seq_interval newsi = gap[ii];
    newsi.SetIntervalId( 0 );
    allw.push_back( newsi );
  }
  sort( allw.begin( ), allw.end( ) );

  // Table with info.
  vec< vec<String> > table;
  table.reserve( allw.size( ) );
  
  vec<seq_interval>::const_iterator it;
  vec< pair<int,int64_t> >::const_iterator it2;
  for (int ii=0; ii<allw.isize( ); ii++) {
    vec<String> tentry( 4 );
    int beg = allw[ii].Begin( );
    int end = allw[ii].End( );
    bool isgap = allw[ii].IntervalId( ) == 0;

    // Readjust to base coordinates (for printing only!)
    if ( isgap && beg > 0 ) beg += ( K - 1 );
    if ( ( isgap && end == super_klen ) || ( ! isgap ) ) end += ( K - 1 );
    
    // tentry[1] and tentry[3] are filled down below if this is a unitig.
    tentry[0] = ToString( end - beg );
    tentry[1] = isgap ? "gap" : "";
    tentry[2] = "[" + ToString( beg ) + ", " + ToString( end ) + ")";
    tentry[3] = "";
    
    // Adjust gap entry if gap overlaps reference gap.
    if ( isgap ) {
      it = lower_bound( super_gaps->begin( ), super_gaps->end( ), allw[ii] );
      int lbpos = distance( super_gaps->begin( ), it );
      bool overlap = false;
      if ( lbpos < super_gaps->isize( ) )
	overlap = allw[ii].HasOverlapWith( (*super_gaps)[lbpos] );
      if ( !overlap && lbpos > 0 )
	overlap = allw[ii].HasOverlapWith( (*super_gaps)[lbpos-1] );
      if ( overlap )
	tentry[1] += " (overlaps ref gap)";
    }

    // Fill tentry[1] and tentry[3] for unitigs.
    else {
      pair<int,int64_t> fake = make_pair( allw[ii].Begin( ), -1 );
      it2 = lower_bound( start2id->begin( ), start2id->end( ), fake );
      int pos = distance( start2id->begin( ), it2 );

      // Default values.
      tentry[1] = "error";
      tentry[3] = "error";

      // Find interval in start2id.
      int s2i_begin = -1;
      int s2i_last = -1;
      for (int jj=pos; jj<start2id->isize( ); jj++) {
	int start = (*start2id)[jj].first;
	if ( start < allw[ii].Begin( ) || start >= allw[ii].End( ) ) break;
	if ( s2i_begin < 0 ) s2i_begin = jj;
	s2i_last = jj;
      }

      // Adjust tentry[1] and [3].
      if ( s2i_begin > -1 ) {
	int vleft = (*to_left)[ (*start2id)[s2i_begin].second ];
	int vright = (*to_right)[ (*start2id)[s2i_last].second ];
	int extleft = ( graphE->ToEdgeObj( vleft ) ).isize( );
	int extright = ( graphE->FromEdgeObj( vright ) ).isize( );

	tentry[1] = "[" + ToString( extleft ) + "]";
	tentry[1] += "[" + ToString( 1 + s2i_last - s2i_begin ) + " edge";
	if ( 1 + s2i_last - s2i_begin > 1 ) tentry[1] += "s";
	tentry[1] += "]";
	tentry[1] += "[" + ToString( extright ) + "]";
	  
	tentry[3] = "v" + ToString( vleft );
	tentry[3] += " " + ToString( (*start2id)[s2i_begin].second );
	if ( 1 + s2i_last - s2i_begin > 2 )
	  tentry[3] += " ...";
	if ( 1 + s2i_last - s2i_begin > 1 )
	  tentry[3] += " " + ToString( (*start2id)[s2i_last].second );
	tentry[3] += " v" + ToString( vright );
      }
      
    }
    
    table.push_back( tentry );
  }

  // Print table.
  out << "\n";
  String justify = "rlll";
  PrintTabular( out, table, 3, justify );
  
}
