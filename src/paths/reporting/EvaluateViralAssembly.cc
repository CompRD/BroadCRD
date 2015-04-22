///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "CoverageAnalyzer.h"
#include "PairsManager.h"
#include "SeqInterval.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"
#include "paths/reporting/CHits.h"
#include "math/Functions.h"
#include "util/RunCommand.h"
// MakeDepend: dependency MakeLookupTable
// MakeDepend: dependency QueryLookupTable
// MakeDepend: dependency RunGetAlignsFast

/**
 * EvaluateViralAssembly
 *
 * Basic assembly stats, geared toward viral assemblies. Alignments
 * and other data will be cached.
 *
 * Input:
 *  <REF_HEAD>.fasta                  (reference fasta)
 *  <RUN_DIR>/reads.{fastb,pairs}     (input reads)
 *  <SUB_DIR>/contigs.fastb           (input contigs)
 *
 * Output: various cached aligns (mostly saved in <SUB_DIR>), plus
 *  info statistics (sent to either cout or file).
 *
 * MIN_CLEN: report long contigs aligning on ref
 * ARCHIVE: send log to <SUB_DIR>/EvaluateViralAssembly.log rather than cout
 * FORCE: do not use cached alignments or data
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( RUN_DIR );
  CommandArgument_String( REF_HEAD );
  CommandArgument_String( SUB_DIR );
  CommandArgument_String_OrDefault( READS, "reads" );
  CommandArgument_String_OrDefault( CONTIGS, "contigs" );
  CommandArgument_Int_OrDefault( MIN_CLEN, 250 );
  CommandArgument_Bool_OrDefault( ARCHIVE, True );
  CommandArgument_Bool_OrDefault( FORCE, False );
  EndCommandArguments;
  
  // File names.
  String sub_dir = RUN_DIR + "/" + SUB_DIR;
  String tmp_dir = sub_dir + "/tmp";

  String readsOnContigs = READS + "_on_" + CONTIGS;
  String readsOnRef = READS + "_on_genome";
  
  String refFasta_file = REF_HEAD + ".fasta";
  String refBases_file = REF_HEAD + ".fastb";
  String refLookup_file = REF_HEAD + ".lookup";
  String readsBases_file = RUN_DIR + "/" + READS + ".fastb";
  String readsPairs_file = RUN_DIR + "/" + READS + ".pairs";
  String contigsHead = RUN_DIR + "/" + SUB_DIR + "/" + CONTIGS;
  String contigsBases_file = contigsHead + ".fastb";
  String contigsLookup_file = contigsHead + ".lookup";
  String contigsOnRef_qlt_file = contigsHead + ".qlt";
  String readsOnRef_qlt_file = sub_dir + "/" + readsOnRef + ".qlt";
  String readsOnContigs_qlt_file = sub_dir + "/" + readsOnContigs + ".qlt";
  String log_file = sub_dir + "/EvaluateViralAssembly.log";

  vec<String> needed;
  needed.push_back( refFasta_file );
  needed.push_back( readsBases_file );
  needed.push_back( readsPairs_file );
  needed.push_back( contigsBases_file );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;
  
  Mkpath( tmp_dir );
  
  // Log stream.
  ofstream log_stream;
  if ( ARCHIVE ) log_stream.open( log_file.c_str( ) );
  ostream &out = ARCHIVE ? log_stream : * (ostream *) &cout;
  if ( ARCHIVE ) {
    PrintCommandPretty( out );
    cout << "Sending log to " << log_file << "\n" << endl;
  }
  
  // Load reference.
  if ( FORCE || ! IsRegularFile( refLookup_file ) ) {
    String theCommand
      = "MakeLookupTable QUIET=True SOURCE=" + refFasta_file
      + " OUT_HEAD=" + REF_HEAD;
    RunCommand( theCommand );
  }
  vecbvec ref( refBases_file );
  vec<int> reflens;
  reflens.reserve( ref.size( ) );
  longlong tot_reflen = 0;
  for (int ii=0; ii<(int)ref.size( ); ii++) {
    reflens.push_back( ref[ii].size( ) );
    tot_reflen += ref[ii].size( );
  }
  
  // Load reads, compute read coverage.
  vecbvec reads( readsBases_file );
  int tot_readlen = 0;
  for (int ii=0; ii<(int)reads.size( ); ii++)
    tot_readlen += reads[ii].size( );
  double cov = SafeQuotient( tot_readlen, tot_reflen );

  out << "READS IN INPUT\n"
      << "\n" << ToStringAddCommas( reads.size( ) ) << " reads in input ("
      << ToString( cov, 1 ) << "x, based on a genome size of\n "
      << ToStringAddCommas( tot_reflen ) << " bases, as computed from the "
      << "given reference).\n" << endl;
  
  // Keep track of which contigs are mapped to reference.
  vecbvec contigs( contigsBases_file );
  vec<bool> cg_select( contigs.size( ), false );
  
  // Aligns of contigs on reference.
  vec<int> cglens;
  cglens.reserve( contigs.size( ) );
  for (int ii=0; ii<(int)contigs.size( ); ii++)
    cglens.push_back( contigs[ii].size( ) );

  // Load aligns of contigs on reference.
  if ( FORCE || ! IsRegularFile( contigsOnRef_qlt_file ) ) {
    String args1 = "K=12 MM=48 MF=500 SH=True MC=0.15 AI=True";
    String args2 = "PARSEABLE=True PRINT_MEMORY_USAGE=True";
    String args3 = "LIST_UNPLACED_BY_PASS=False";
    
    String theCommand
      = "QueryLookupTable " + args1 + " " + args2 + " " + args3
      + " L=" + refLookup_file
      + " SEQS=" + contigsBases_file
      + " OUTPUT=" + contigsOnRef_qlt_file
      + " TMP_DIR=" + tmp_dir;
    RunCommand( theCommand );
  }
  
  vec<look_align> hits;
  LoadLookAligns( contigsOnRef_qlt_file, hits );
  order_lookalign_TargetBegin sorter;
  sort( hits.begin( ), hits.end( ), sorter );

  // Big contigs (skipping subsumed chunks).
  vec< pair<int,int> > len2hit;
  len2hit.reserve( hits.size( ) );
  for (int ii=0; ii<hits.isize( ); ii++) {
    const look_align &al = hits[ii];
    if ( (int)al.query_length < MIN_CLEN ) continue;
    int len = al.a.Pos2( ) - al.a.pos2( );
    len2hit.push_back( make_pair( len, ii ) );
  }
  sort( len2hit.rbegin( ), len2hit.rend( ) );
  vec<seq_interval> big_contigs;
  for (int ii=0; ii<len2hit.isize( ); ii++) {
    const int &hit_id = len2hit[ii].second;
    const look_align &al = hits[hit_id];
    int seq_id = al.target_id;
    int begin = al.a.pos2( );
    int end = al.a.Pos2( );
    bool embedded = false;
    for (int jj=0; jj<big_contigs.isize( ); jj++) {
      if ( big_contigs[jj].SeqId( ) != seq_id ) continue;
      if ( big_contigs[jj].Begin( ) > begin ) continue;
      if ( big_contigs[jj].End( ) < end ) continue;
      embedded = true;
      break;
    }
    if ( ! embedded )
      big_contigs.push_back( seq_interval( hit_id, seq_id, begin, end ) );
  }
  vec<int> big_lens;
  vec<bool> is_big( hits.size( ), false );
  for (int ii=0; ii<big_contigs.isize( ); ii++) {
    is_big[ big_contigs[ii].IntervalId( ) ] = true;
    big_lens.push_back( hits[ big_contigs[ii].IntervalId( ) ].query_length );
  }
  sort( big_lens.begin( ), big_lens.end( ) );
  int big_n50 = N50( big_lens );

  // Table of contigs >= MIN_CLEN aligned to ref.
  vec< vec<String> > table;
  vec<String> line = MkVec( String( "cg_id" ),
			    String( "align_on_contig" ),
			    String( "reference_id" ),
			    String( "align_on_reference" ),
			    String( "big" ) );
  table.push_back( line );
  for (int ii=0; ii<hits.isize( ); ii++) {
    const look_align &al = hits[ii];
    if ( (int)al.query_length < MIN_CLEN ) continue;
    cg_select[al.query_id] = true;
    String al_cg
      = "[" + ToString( al.a.pos1( ) )
      + ", " + ToString( al.a.Pos1( ) )
      + ")_" + ToString( al.query_length );
    String al_ref
      = "[" + ToString( al.a.pos2( ) )
      + ", " + ToString( al.a.Pos2( ) )
      + ")_" + ToString( al.target_length );
    String t_id = " on t" + ToString( al.target_id );
    line.clear( );
    line.push_back( ToString( al.query_id ),
		    al_cg,
		    String( al.Rc1( ) ? "-" : "+" ) + t_id,
		    al_ref,
		    ( is_big[ii] ? "#####" : "" ) );
    table.push_back( line );
  }

  out << "ALIGNS OF CONTIGS >= " << MIN_CLEN << " BASES ON REFERENCE\n\n";
  PrintTabular( out, table, 2, "rrrrl" );
  out << endl;
  
  // Load pairs (at this time only supported one library).
  PairsManager pairs( readsPairs_file );
  ForceAssertEq( (int)pairs.nLibraries( ), 1 );

  // Generate lookup table of contigs.
  if ( FORCE || ! IsRegularFile( contigsLookup_file ) ) {
    String theCommand
      = "MakeLookupTable LOOKUP_ONLY=True QUIET=True SOURCE="+contigsBases_file
      + " OUT_HEAD=" + contigsHead;
    RunCommand( theCommand );
  }
  
  // Align reads on contigs.
  if ( FORCE || ! IsRegularFile( readsOnContigs_qlt_file ) ) {
    String theCommand
      = "RunGetAlignsFast K=20 FASTB_FILE=" + readsBases_file
      + " LOOKUP_FILE=" + contigsLookup_file
      + " ALIGNS_FILE=" + readsOnContigs_qlt_file
      + " USE_CACHE=" + ( FORCE ? "False" : "True" )
      + " TMP_DIR=" + tmp_dir;
    RunCommand( theCommand );
  }

  // Report regions on assembly with no physical coverage.
  vec<look_align> roc_hits;
  LoadLookAligns( readsOnContigs_qlt_file, roc_hits );
    
  vec<seq_interval> nocov_contigs;
  {
    out << "ESTIMATED LIBRARY STATS (ON ASSEMBLY)\n\n";
    CHits analyzer( roc_hits, pairs, cglens, &out );
    
    int lib_id = 0;
    int min_cov = 1;
    analyzer.CoveragesAtMost( lib_id, min_cov - 1, nocov_contigs );
    
    vec< vec<String> > table;
    vec<String> line = MkVec( String( "cg_id" ),
			      String( "beg" ),
			      String( "end" ),
			      String( "cg_len" ) );
    table.push_back( line );
    for (int ii=0; ii<nocov_contigs.isize( ); ii++) {
      if ( ! cg_select[ nocov_contigs[ii].SeqId( ) ] ) continue;
      line.clear( );
      line = MkVec( ToString( nocov_contigs[ii].SeqId( ) ),
		    ToString( nocov_contigs[ii].Begin( ) ),
		    ToString( nocov_contigs[ii].End( ) ),
		    ToString( cglens[ nocov_contigs[ii].SeqId( ) ] ) );
      table.push_back( line );
    }
    out << "INTERVALS WITH NO PHYSICAL COVERAGE (CONTIGS)" << "\n\n";
    PrintTabular( out, table, 2, "rrrr" );
    out << endl;
  }
  
  // Align reads on reference.
  if ( FORCE || ! IsRegularFile( readsOnRef_qlt_file ) ) {
    String theCommand
      = "RunGetAlignsFast K=20 FASTB_FILE=" + readsBases_file
      + " LOOKUP_FILE=" + refLookup_file
      + " ALIGNS_FILE=" + readsOnRef_qlt_file
      + " USE_CACHE=" + ( FORCE ? "False" : "True" )
      + " TMP_DIR=" + RUN_DIR;
    RunCommand( theCommand );
  }

  // Report regions on reference with no physical coverage.
  vec<look_align> rog_hits;
  LoadLookAligns( readsOnRef_qlt_file, rog_hits );
  
  vec<seq_interval> cov_ref;
  vec<seq_interval> nocov_ref;
  {
    out << "ESTIMATED LIBRARY STATS (ON REFERENCE)\n\n";
    CHits analyzer( rog_hits, pairs, reflens, &out );
  
    int lib_id = 0;
    int min_cov = 1;
    analyzer.CoveragesAtMost( lib_id, min_cov - 1, nocov_ref );
    analyzer.CoveragesAtLeast( lib_id, min_cov, cov_ref );
    
    vec< vec<String> > table;
    vec<String> line = MkVec( String( "ref_id" ),
			      String( "beg" ),
			      String( "end" ) );
    table.push_back( line );
    for (int ii=0; ii<nocov_ref.isize( ); ii++) {
      line.clear( );
      line = MkVec( ToString( nocov_ref[ii].SeqId( ) ),
		    ToString( nocov_ref[ii].Begin( ) ),
		    ToString( nocov_ref[ii].End( ) ) );
      table.push_back( line );
    }
    out << "INTERVALS WITH NO PHYSICAL COVERAGE (REFERENCE)" << "\n\n";
    PrintTabular( out, table, 2, "rrr" );
    out << endl;
  }
  
  // Parse selected aligns for: %covered, %error rate, misassemblies.
  int n_misassemblies = 0;
  int n_errors = 0;
  int aligned_len = 0;
  vec<seq_interval> covered_regions;
  for (int ii=0; ii<hits.isize( ); ii++) {
    const look_align &al = hits[ii];
    if ( ! cg_select[al.query_id] ) continue;
    
    // Covered portion.
    seq_interval creg( ii, al.target_id, al.a.pos2( ), al.a.Pos2( ) );
    covered_regions.push_back( creg );

    // Error rate.
    aligned_len += al.a.Pos1( ) - al.a.pos1( );
    n_errors += al.Errors( );

    // Misassemblies.
    int clen = al.query_length;
    int tlen = al.target_length;
    if ( al.a.pos1( ) != 0 && al.a.pos2( ) != 0 ) n_misassemblies++;
    if ( al.a.Pos1( ) != clen && al.a.Pos2( ) != tlen ) n_misassemblies++;
  }

  // Regions with read coverage.
  CoverageAnalyzer canalyzer( covered_regions, &reflens );
  vec<seq_interval> read_andcg_cov;
  canalyzer.GetCoveragesAtLeast( 1, read_andcg_cov );
  
  // Simple trick to find regions on ref covered by reads and contigs.
  vec<seq_interval> temp_cov;
  temp_cov.reserve( 2 * read_andcg_cov.size( ) );
  for (int ii=0; ii<read_andcg_cov.isize( ); ii++) {
    temp_cov.push_back( read_andcg_cov[ii] );
    temp_cov.push_back( read_andcg_cov[ii] );
  }
  swap( temp_cov, read_andcg_cov );
  copy( cov_ref.begin( ), cov_ref.end( ), back_inserter( read_andcg_cov ) );
  
  vec<seq_interval> temp_cov2;
  vec<seq_interval> read_andcg_nocov;
  CoverageAnalyzer canalyzer2( read_andcg_cov, &reflens );
  canalyzer2.GetCoveragesExactly( 3 , temp_cov2 );
  canalyzer2.GetCoveragesExactly( 1 , read_andcg_nocov );
  swap( temp_cov2, read_andcg_cov );
  int read_andcg_len = 0;
  for (int ii=0; ii<read_andcg_cov.isize( ); ii++)
    read_andcg_len += read_andcg_cov[ii].Length( );

  // Regions with read coverage but no contig coverage.
  out << "REGIONS ON REFERENCE WITH READ COVERAGE BUT NO CONTIG COVERAGE\n\n";
  for (int ii=0; ii<read_andcg_nocov.isize( ); ii++)
    out << "t" << read_andcg_nocov[ii].SeqId( )
	<< "\t[" << read_andcg_nocov[ii].Begin( )
	<< ", " << read_andcg_nocov[ii].End( )
	<< ")\n";
  out << endl;

  int cov_ref_len = 0;
  for (int ii=0; ii<cov_ref.isize( ); ii++)
    cov_ref_len += cov_ref[ii].Length( );
  double cov_ref_ratio = 100.0 * SafeQuotient( cov_ref_len, tot_reflen );
  double read_andcg_ratio = 100.0 * SafeQuotient( read_andcg_len, cov_ref_len );
  
  // Print summary.
  out << "ERRORE RATE, COVERAGE, AND MISASSEMBLIES (BY CONTIGS >= "
      << MIN_CLEN << " BASES)\n\n"
      << ToString( cov_ref_ratio, 1 ) << "%   -   "
      << "percent of reference covered by reads\n"
      << ToString( read_andcg_ratio, 1 ) << "%   -   "
      << "percent of reference covered by reads, that is also\n"
      << "  covered by contigs\n"
      << ToString( 100. * SafeQuotient(n_errors,aligned_len), 1 ) << "%   -   "
      << "percent error rate (computed as 100.0 times the total number\n"
      << "  of errors, divided by the total aligned length of contigs)\n"
      << n_misassemblies << "   -   "
      << "number of misassemblies (computed as the number of times a contig\n"
      << "  aligns the reference with an hanging end).\n"
      << endl;
  
  // A csv list (can be grepped and inserted in a spreadsheet).
  out << "CSV SUMMARY TABLE\n\n"
      << "Legend\n"
      << "  n_reads:   number of reads in input\n"
      << "  r_cov:     read coverage (based on a genome size of "
      << ToString( tot_reflen ) << " bases, as\n"
      << "             computed from the given reference)\n"
      << "  n_contigs: number of big contigs (ie contigs >= " << MIN_CLEN
      << " bases,\n"
      << "             after removing subsumed contigs)\n"
      << "  N50:       N50 of big contigs\n"
      << "  rd_cov:    \"reads-on-reference coverage\" (percent of reference\n"
      << "             covered by reads)\n"
      << "  cg_cov:    \"contigs-on-reference coverage\" (percent of"
      << " reference\n"
      << "             covered by reads, that is also covered by contigs)\n"
      << "  err_rate:  percent error rate (computed as 100.0 times the total\n"
      << "             number of errors, divided by the total aligned length\n"
      << "             of contigs)\n"
      << "  n_mis:     number of misassemblies (computed as the number of\n"
      << "             times a contig aligns the reference with an hanging\n"
      << "             end)\n"
      << "  run_dir:   RUN_DIR arg\n"
      << "\n"
      << "csv_labels,"
      << "n_reads,"
      << "r_cov (x),"
      << "n_contigs,"
      << "N50,"
      << "rd_cov (%),"
      << "cg_cov (%),"
      << "err_rate (%),"
      << "n_mis,"
      << "run_dir\n"
      << "csv_data,"
      << reads.size( ) << ","
      << ToString( cov, 1 ) << ","
      << big_contigs.size( ) << ","
      << big_n50 << ","
      << ToString( cov_ref_ratio, 1 ) << ","
      << ToString( read_andcg_ratio, 1 ) << ","
      << ToString( 100. * SafeQuotient( n_errors, aligned_len ), 1 ) << ","
      << n_misassemblies << ","
      << RUN_DIR << "\n"
      << endl;
  
  // Done.
  cout << Date( ) << ": EvaluateViralAssembly done" << endl;
  
}

