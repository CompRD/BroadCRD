///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "CoverageAnalyzer.h"
#include "Superb.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "math/NStatsTools.h"
#include "util/RunCommand.h"

void AddLine( const String &name,
	      vec<int> &lens,
	      vec< vec<String> > &table,
	      ostream *blog = 0 );

void EvaluateAligns( const String &aligns_file,
		     const String &ref_fastb_file,
		     ostream &log,
		     ostream &blog );

/**
 * BuildAssistedReport
 *
 * Basic report utility for assisted assemblies. It optionally
 * evaluates the assembly against a reference, if REF_HEAD is not
 * empty.
 *
 * NOTE (if REF_HEAD is given): you must run EvalScaffolds to generate
 *   the aligns of the contigs onto the reference.
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( ASSEMBLY );
  CommandArgument_String_OrDefault( REF_HEAD, "" );
  CommandArgument_String_OrDefault( OUT_BASE, ASSEMBLY + ".report" );
  EndCommandArguments;

  // File names.
  String supers_file = ASSEMBLY + ".superb";
  String chains_file = ASSEMBLY + ".Eval/chains.qlt";
  String ref_fastb_file = REF_HEAD + ".fastb";

  String log_file = OUT_BASE;
  String brief_file = OUT_BASE + ".brief";

  vec<String> needed;
  needed.push_back( supers_file );
  if ( REF_HEAD != "" ) needed.push_back( ref_fastb_file );
  if ( REF_HEAD != "" ) needed.push_back( chains_file );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;

  // Output stream.
  ofstream log( log_file.c_str( ) );
  ofstream blog( brief_file.c_str( ) );
  PrintCommandPretty( log );
  blog << "ASSEMBLY\t" << ASSEMBLY << "\n";
  
  // Load supers.
  vec<superb> supers;
  ReadSuperbs( supers_file, supers );

  // Lengths containers.
  vec<int> slens;          // super lengths
  vec<int> clens;          // contig lenghts
  vec<int> neg_glens;      // gap lengths (<0)
  vec<int> nonneg_glens;   // gap lengths (>=0)
  vec<int> all_glens;      // gap lengths (all)

  // Various sizes, for reserving.
  int n_supers = supers.size( );
  int n_contigs = 0;
  int n_neg_gaps = 0;
  int n_nonneg_gaps = 0;
  for (int ii=0; ii<n_supers; ii++) {
    int nc = supers[ii].Ntigs( );
    n_contigs += nc;
    for (int jj=0; jj<nc-1; jj++) {
      int gap = supers[ii].Gap( jj );
      if ( gap < 0 ) n_neg_gaps++;
      else n_nonneg_gaps++;
    }
  }
  
  // Reserve.
  slens.reserve( n_supers );
  clens.reserve( n_contigs );
  neg_glens.reserve( n_neg_gaps );
  nonneg_glens.reserve( n_nonneg_gaps );
  all_glens.reserve( n_neg_gaps + n_nonneg_gaps );

  // Fill lengths containers (they will be sorted in AddLine below).
  for (int ii=0; ii<n_supers; ii++) {
    slens.push_back( supers[ii].TrueLength( ) );
    int nc = supers[ii].Ntigs( );
    for (int jj=0; jj<nc; jj++) {
      clens.push_back( supers[ii].Len( jj ) );
      if ( jj == nc - 1 ) continue;
      int gap = supers[ii].Gap( jj );
      all_glens.push_back( gap );
      if ( gap < 0 ) neg_glens.push_back( -gap );
      else nonneg_glens.push_back( gap );
    }
  }
  
  // Build and print core assembly stats.
  vec< vec<String> > table;
  
  table.push_back( MkVec( String( "" ),
			  String( "#" ),
			  String( "range" ),
			  String( "median" ),
			  String( "N50" ),
			  String( "total length" ) ) );

  AddLine( "scaffolds", slens, table, &blog );
  AddLine( "contigs", clens, table, &blog );
  AddLine( "gaps (<0)", neg_glens, table );
  AddLine( "gaps (>=0)", nonneg_glens, table );
  AddLine( "gaps (all)", all_glens, table );
  
  log << "ASSEMBLY STATS\n\n";
  PrintTabular( log, table, 3, "lrrrrr" );
  log << endl;

  // Evaluate aligns on reference.
  if ( REF_HEAD != "" )
    EvaluateAligns( chains_file, ref_fastb_file, log, blog );

  // Done.
  log << Date( ) << ": BuildAssistedReport done" << endl;
  log.close( );
  blog.close( );

}



////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//     FUNCTIONS                                                              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////



/**
 * AddLine
 *
 * Add a line to the assembly stats table (it will sort lens).
 */
void AddLine( const String &name,
	      vec<int> &lens,
	      vec< vec<String> > &table,
	      ostream *blog )
{
  sort( lens.begin( ), lens.end( ) );

  String s_count = ToString( lens.size( ) );
  
  String s_range = "na";
  if ( lens.size( ) > 0 ) {
    String strF = ToString( lens[0] );
    String strL = ToString( lens[ lens.size( )-1 ] );
    s_range = "[" + strF + ", " + strL + "]";
  }

  String s_med = "na";
  if ( lens.size( ) > 0 )
    s_med = ToString( Median( lens ) );

  String s_n50 = "na";
  if ( lens.size( ) > 0 && lens[0] >= 0 ) {
    vec<int> clean_lens;
    clean_lens.reserve( lens.size( ) );
    for (int ii=0; ii<lens.isize( ); ii++)
      if ( lens[ii] > 0 ) clean_lens.push_back( lens[ii] );
    if ( clean_lens.size( ) > 0 )
      s_n50 = ToString( N50( clean_lens ) );
  }
  
  String s_tot = "na";
  if( lens.size( ) > 0 )
    s_tot = ToString( BigSum( lens ) );
  
  table.push_back( MkVec( name, s_count, s_range, s_med, s_n50, s_tot ) );
  
  if ( blog ) {
    *blog << "N" << name << "\t" << s_count << "\n"
	  << "N50" << name << "\t" << s_n50 << "\n"
	  << "TOTLEN" << name << "\t" << s_tot << "\n";
  }
  
}

/**
 * EvaluateAligns
 *
 * Load chains of aligns of contigs on reference, and generate overall
 * statistics.
 */
void EvaluateAligns( const String &aligns_file,
		     const String &ref_fastb_file,
		     ostream &log,
		     ostream &blog )
{
  // Load data.
  vecbvec ref( ref_fastb_file );
  
  vec<look_align> aligns;
  LoadLookAligns( aligns_file, aligns );
  
  // Reference lengths.
  vec<int> rlens;
  rlens.reserve( ref.size( ) );
  for (size_t ii=0; ii<ref.size( ); ii++)
    rlens.push_back( ref[ii].size( ) );
  
  longlong tot_rlen = BigSum( rlens );
  
  // Analyze coverage.
  longlong tot_mutations = 0;
  longlong tot_indels = 0;
  vec<seq_interval> raw_sis;
  raw_sis.reserve( aligns.size( ) );
  for (size_t ii=0; ii<aligns.size( ); ii++) {
    int sid = aligns[ii].target_id;
    int beg = aligns[ii].pos2( );
    int end = aligns[ii].Pos2( );
    raw_sis.push_back( seq_interval( (int)ii, sid, beg, end ) );
    tot_mutations += aligns[ii].Mutations( );
    tot_indels += aligns[ii].Indels( );
  }
  
  CoverageAnalyzer cov( raw_sis, &rlens );

  vec<seq_interval> uncovered_sis;
  cov.GetCoveragesExactly( 0, uncovered_sis );
  longlong tot_uncovered_len = 0;
  for (size_t ii=0; ii<uncovered_sis.size( ); ii++)
    tot_uncovered_len += uncovered_sis[ii].Length( );

  vec<seq_interval> covered_sis;
  cov.GetCoveragesAtLeast( 1, covered_sis );
  longlong tot_covered_len = 0;
  for (size_t ii=0; ii<covered_sis.size( ); ii++)
    tot_covered_len += covered_sis[ii].Length( );
  
  ForceAssertEq( tot_covered_len + tot_uncovered_len, tot_rlen );
  longlong tot_errors = tot_mutations + tot_indels;
  double ratio_cov = SafeQuotient( tot_covered_len, tot_rlen );
  double ratio_mut = SafeQuotient( tot_mutations, tot_covered_len );
  double ratio_indels = SafeQuotient( tot_indels, tot_covered_len );
  double ratio_err = SafeQuotient( tot_errors, tot_covered_len );

  // Coverage table.
  vec< vec<String> > table;

  table.push_back( MkVec( String( "" ),
			  String( "bases" ),
			  String( "% of genome" ) ) );
  
  table.push_back( MkVec( String( "coverage" ),
			  ToStringAddCommas( tot_covered_len ),
			  ToString( 100. * ratio_cov, 2 ) ) );

  log << "COVERAGE OF REFERENCE (based on a genome size of "
      << ToStringAddCommas( tot_rlen ) << " bases)\n\n";
  PrintTabular( log, table, 3, "lrr" );
  log << endl;

  // Errors table.
  table.clear( );

  table.push_back( MkVec( String( "" ),
			  String( "bases" ),
			  String( "% of covered portion" ) ) );
  
  table.push_back( MkVec( String( "errors (mutations only)" ),
			  ToStringAddCommas( tot_mutations ),
			  ToString( 100. * ratio_mut, 2 ) ) );
  
  table.push_back( MkVec( String( "errors (indels only)" ),
			  ToStringAddCommas( tot_indels ),
			  ToString( 100. * ratio_indels, 2 ) ) );

  table.push_back( MkVec( String( "erors (total)" ),
			  ToStringAddCommas( tot_errors ),
			  ToString( 100. * ratio_err, 2 ) ) );

  log << "ERRORS STATISTICS (based on a covered portion of "
      << ToStringAddCommas( tot_covered_len ) << " bases)\n\n";
  PrintTabular( log, table, 3, "lrr" );
  log << endl;
  
  // Generate brief output file.
  blog << "COVERED\t" << tot_covered_len << "\n"
       << "MUTATIONS\t" << tot_mutations << "\n"
       << "INDELS\t" << tot_indels << "\n"
       << endl;
  
}
