///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


/* EvalHyperLG: A HyperKmerPath evaluation tool.
 * Designed for use with the RunAllPathsLG pipeline.
 *
 * Adapted from the original EvalHyper via massive disentanglement and cleanup.
 *
 * Josh Burton
 * June 2009
 *
 ******************************************************************************/

#include "Basevector.h"
#include "Bitvector.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "lookup/LookAlign.h"
#include "paths/AlignHyperKmerPath.h"
#include "paths/EvalUtilsLG.h" // main utils module
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"












int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  
  CommandDoc( "Evaluates the HyperKmerPath at <PRE>/<DATA>/<RUN>/ASSEMBLIES/<SUBDIR>/<HYPER>.  Set the EVAL flags to control what evaluation is performed.  Each flag EVAL_FOO creates a file <SUBDIR>/EvalHyperLG/EvalHyperLG.FOO.out." );
  
  // Directory and file names.
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String( SUBDIR );
  
  CommandArgument_String_OrDefault( READS, "all_reads" );
  CommandArgument_String_OrDefault( WRUN, "" );
  CommandArgument_String_OrDefault( HYPER, "hyper" );
  CommandArgument_String_OrDefault_Doc( REFTIGS, "all",
      "Choose reference supercontigs for evaluation i.e. \"all\" or \"1,5,18\"" );
  CommandArgument_Bool_OrDefault_Doc( FILTER_ALIGNS, True, "Filter alignments and write file hyper.trusted_paths" );
  
  // EVAL flags.
  // Each of these flags, when set to True, causes a file to be created in
  // sub_dir/EvalHyperLG, with a name implied by the flag's name.  For example,
  // if EVAL_LIB_STATS=True, a file is written:
  // <sub_dir>/EvalHyperLG/EvalHyperLG.lib_stats.out
  // Lastly, these files are concatenated into EvalHyperLG/EvalHyperLG.all.out.
  CommandArgument_Bool_OrDefault_Doc( EVAL_SUMMARY, True, "Basic data: reference coverage, graph statistics, errors, gaps, etc.  Also writes to stdout." );
  CommandArgument_Bool_OrDefault_Doc( EVAL_LIB_STATS, True, "Library insert-size statistics." );
  CommandArgument_Bool_OrDefault_Doc( EVAL_ALIGN_TO_REF, True, "Alignment of HyperKmerPath to reference." );
  CommandArgument_String_OrDefault_Doc( LINEAR_SCAFFOLDS, "", "Evaluate scaffolds using ScaffoldAccuracy module using this file" );
  
  EndCommandArguments;
  
  
  
  // Set up directories.
  const String data_dir = PRE + "/" + DATA;
  const String run_dir = PRE + "/" + DATA + "/" + RUN;
  const String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  const String wrun_dir = sub_dir + "/" + WRUN;
  String genome_head = "genome";
  
  String output_dir = sub_dir + "/EvalHyperLG";
  String hyper_suffix = HYPER.After( "hyper" );
  if ( hyper_suffix.StartsWith( "_" ) ) hyper_suffix = hyper_suffix.After( "_" );
  if ( hyper_suffix != "" ) output_dir += hyper_suffix;
  Mkdir777( output_dir );
  String hyper_file = sub_dir + "/" + HYPER;
  
  // Determine ploidy.
  const int ploidy = FirstLineOfFile( data_dir + "/ploidy" ).Int( );
  const Bool diploid = ( ploidy > 1 );
  
  
  /*****************************************************************************
   *
   *                          LOAD DATA
   *
   *   Some data structures are only needed if certain EVAL flags have been set.
   *   We save time and memory by determining which files are necessary, and
   *   loading only those.
   *
   *   Data structures that are declared here but not loaded from file will
   *   remain as dummy variables, which is okay.
   *
   ****************************************************************************/
  
  // Data structure declarations.
  vecbasevector genome, genome_diploid;
  vecbitvector genome_amb;
  PairsManager pairs;
  HyperKmerPath hkp;
  KmerBaseBroker * kbb = NULL;
  vec<look_align> aligns;
  vec< vec<int> > aligns_index;
  vec<TrustedPath> trusted_paths;
  
  
  cout << Date( ) << ": Loading files and creating data structures..." << endl;

  vec<int> ref_scontigs;
  if ( REFTIGS != "all" ){
    vec<char> separators; 
    separators.push_back( ',');
    vec<String> sref_scontigs;
    Tokenize( REFTIGS, separators, sref_scontigs );
    for ( unsigned i = 0; i < sref_scontigs.size(); i++ ){
      ref_scontigs.push_back( atoi( sref_scontigs[i].c_str() ) );
    }
    cout << "Reference supercontigs for use:\n ";
    ref_scontigs.Print( cout ); cout << endl;
    
    vecbasevector subgenome;
    vecbitvector subgenome_amb;
    subgenome.Read( data_dir + "/" + genome_head + ".fastb", ref_scontigs );
    if ( IsRegularFile( data_dir + "/" + genome_head + ".fastamb" ) )
      subgenome_amb.Read( data_dir + "/" + genome_head + ".fastamb", ref_scontigs );

    subgenome.WriteAll( data_dir + "/subgenome.fastb" );
    if ( IsRegularFile( data_dir + "/" + genome_head + ".fastamb" ) )
      subgenome_amb.WriteAll( data_dir + "/subgenome.fastamb" );

    SystemSucceed( "MakeLookupTable" + ARG(SOURCE, data_dir + "/subgenome.fastb" )
		   + ARG(OUT_HEAD, data_dir + "/subgenome") + ARG(LOOKUP_ONLY, True)
		   + ARG( NH, True ) );

    genome_head = "subgenome";
  }
  
  if ( EVAL_LIB_STATS ) {
    longlong n_reads = MastervecFileObjectCount( run_dir + "/" + READS + ".fastb" );
    pairs.Read( run_dir + "/" + READS + ".pairs" );
  }
  
  if ( EVAL_SUMMARY || EVAL_ALIGN_TO_REF )
      genome = vecbasevector( data_dir + "/" + genome_head + ".fastb" );
  
  
  if ( EVAL_SUMMARY ) {
    if ( IsRegularFile( data_dir + "/genome.diploid.fastb" ) )
      genome_diploid.ReadAll( data_dir + "/genome.diploid.fastb" );
    if ( IsRegularFile( data_dir + "/" + genome_head + ".fastamb" ) )
	genome_amb.ReadAll( data_dir + "/" + genome_head + ".fastamb" );
      
    
  }
  
  if ( EVAL_SUMMARY || EVAL_ALIGN_TO_REF ) {
    hkp = HyperKmerPath( hyper_file );
    int K = hkp.K( );
    kbb = new KmerBaseBroker( wrun_dir, K );
    
    // Align the HyperKmerPath to reference.  This is time-consuming.
    cout << Date( ) << ": Aligning HyperKmerPath to reference..." << endl;
    AlignHyperKmerPath( hkp, kbb, data_dir + "/" + genome_head, wrun_dir, aligns, aligns_index );
    
    if ( FILTER_ALIGNS ) {
      // Filter alignments.
      cout << Date( ) << ": Filtering alignments..." << endl;
      FilterAligns( hkp, aligns, aligns_index, trusted_paths, -1 );
    }
    
    // Record trusted paths.  If FILTER_ALIGNS=False, trusted_paths will
    // be an empty vector.
    String trusted_paths_file = hyper_file + ".trusted_paths";
    BinaryWriter::writeFile( trusted_paths_file, trusted_paths );
  }
  
  
  
  
  /*****************************************************************************
   *
   *                          RUN THE EVALUATIONS!
   *
   ****************************************************************************/
  
  // Create summary evaluation.  Write it to a file, as well as to stdout.
  // Creates file: EvalHyperLG.summary.out

  if ( EVAL_SUMMARY ) {
    EvalSummary( "stdout",
		 hkp, genome, genome_diploid, genome_amb, aligns, aligns_index );
    EvalSummary( output_dir + "/EvalHyperLG.summary.out",
		 hkp, genome, genome_diploid, genome_amb, aligns, aligns_index );
  }
  
  // Generate library stats. (Creates file: EvalHyperLG.lib_stats.out)
  if ( EVAL_LIB_STATS )
    EvalLibraryStats( output_dir + "/EvalHyperLG.lib_stats.out", pairs );
  
  
  // Write out the alignment of the HyperKmerPath to reference.
  // Creates file: EvalHyperLG.align_to_ref.out

  if ( EVAL_ALIGN_TO_REF ) {
    String outfile = output_dir + "/EvalHyperLG.align_to_ref.out";
    cout << Date( ) << ": EVAL - align to ref" << endl;
    ofstream out( outfile.c_str( ) );
    PrintAlignedHyperKmerPath( out, hkp, kbb, genome, aligns, aligns_index,
			       True, &trusted_paths, False, diploid );
    out.close( );
  }
  
  // Generate scaffold stats. (Creates file: EvalHyperLG.scaffold_accuracy.out)

  if ( LINEAR_SCAFFOLDS != "")
    EvalScaffoldAccuracy( output_dir + "/EvalHyperLG.scaffold_accuracy.out",
			  data_dir + "/../" + genome_head + ".fasta",
			  sub_dir + "/" + LINEAR_SCAFFOLDS );
  
  // Concatenate output files into a single output file.

  String all_output_file = sub_dir + "/report.brief";
  if ( hyper_suffix != "" ) all_output_file += "." + hyper_suffix;
  Remove( all_output_file );
  if ( EVAL_ALIGN_TO_REF )
    System( "cat " + output_dir + "/EvalHyperLG.align_to_ref.out >> "
	    + all_output_file );
  if ( EVAL_SUMMARY )
    System( "cat " + output_dir + "/EvalHyperLG.summary.out >> "
	    + all_output_file );
  if ( LINEAR_SCAFFOLDS != "")
    System( "cat " + output_dir + "/EvalHyperLG.scaffold_accuracy.out >> "
	    + all_output_file );
  if ( EVAL_LIB_STATS )
    System( "cat " + output_dir + "/EvalHyperLG.lib_stats.out > "
            + sub_dir + "/report.data" );
  
  
  // Clean up memory.
  if ( kbb ) delete kbb;

  // Done!  Tell the user where to find output files.

  cout << Date( ) << ": Done!" << endl << endl;
  if ( EVAL_SUMMARY )
    cout << "This report, along with much more information, has been written to files at" << endl;
  else
    cout << "Report files are at" << endl;
  cout << RealPath( output_dir ) << "/EvalHyperLG.*.out." << endl;
}
