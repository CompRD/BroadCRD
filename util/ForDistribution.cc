/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "util/RunCommand.h"

/**
 * ForDistribution
 *
 * Prepares an assembly for distribution. It creates a subdir in DATA,
 * with human-readable files for the assembly specified in SUBDIR.
 *
 * Remark: no data in SUBDIR will be changed (this module does not
 * change the assembly). You may want to run util/FinalizeAssembly
 * before ForDistribution to prepare the assembly for distribution.
 *
 * BASE_OUT: output directory name (relative to DATA)
 * KNOWN_CONTIGS: known fasta (as a file in DATA)
 * SUPER_TO_FASTA: run SuperToFasta
 * READ_TABLES: run CalculateReadStats
 * ASSEMBLY_TABLES: run CalculateMapStats
 * CHROMOSOME_DATA: it generates fasta/qual based on an agp file (supers based)
 * MARKUP_ONLY: only upgrade markup info
 * FORCE: overwrite BASE_OUT if it already exists
 */
int main( int argc, char *argv[] )
{ 
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String( SUBDIR );
  CommandArgument_String_OrDefault( BASE_OUT, "ForDistribution" );
  CommandArgument_String_OrDefault( KNOWN_CONTIGS, "contigs.fasta" );
  CommandArgument_String_OrDefault( REPORT_DIR, "Report" );
  CommandArgument_String_OrDefault( REPORT_PRETTY_NAME, "assembly.ps" );
  CommandArgument_String_OrDefault( AGPFILE, "assembly.agp" );
  CommandArgument_Bool_OrDefault( SUPER_TO_FASTA, True );
  CommandArgument_Bool_OrDefault( READ_TABLES, False );
  CommandArgument_Bool_OrDefault( ASSEMBLY_TABLES, False );
  CommandArgument_Bool_OrDefault( CHROMOSOME_DATA, True );
  CommandArgument_Bool_OrDefault( MARKUP_ONLY, False );
  CommandArgument_Bool_OrDefault( MARKUP, True );
  CommandArgument_Bool_OrDefault( FORCE, False );
  EndCommandArguments;

  // Dir names.
  String data_dir = PRE + "/" + DATA;
  String run_dir = data_dir + "/" + RUN;
  String sub_dir = run_dir + "/" + SUBDIR;
  String distrib_dir = data_dir + "/" + BASE_OUT;

  // File names.
  String log_file = sub_dir + "/ForDistribution.log";
  String sub_agp_file = sub_dir + "/" + AGPFILE;
  String qual_qc_file = sub_dir + "/qual.qc.xml";
  String lowcov_qc_file = sub_dir + "/low_coverage.qc.xml";

  String source_file = distrib_dir + "/source";
  String comm_file = distrib_dir + "/ForDistribution.command";
  String markup_file = distrib_dir + "/assembly.markup";
  String agp_file = distrib_dir + "/assembly.agp";
  String readme_file = distrib_dir + "/README";

  String oneliner_file = "BasicAssemblyOneLiner.out";
  String bastats_file = "BasicAssemblyStats.out";
  String libstats_file = "LibStatsOverview.out";
  String physcov_file = "PhysicalCoverageByLib.out";

  // Early exit if needed.
  if ( IsDirectory( distrib_dir ) && !FORCE ) {
    cout << "Warning! A directory with the name of\n " << distrib_dir
	 << "\nalready exists. You can either use another name,\n"
	 << "or call ForDistribution with FORCE=True\n"
	 << "to overwrite the existing one. Now exit.\n\n";
    return -1;
  }
  
  Mkdir777( distrib_dir );

  // README file.
  ofstream orme( readme_file.c_str( ) );
  orme << "See http://www.broad.mit.edu/wga/arachnewiki/index.php/Output\n"
       << "for details about the files in this directory." << endl;
  orme.close( );

  // Main log file.
  ofstream log( log_file.c_str( ), ios::app );
  PrintCommandPretty( log );
  
  // Record calling command, and source SUBDIR.
  if ( !MARKUP_ONLY ) {
    ofstream out_comm( comm_file.c_str( ) );
    PrintCommandPretty( out_comm );
    out_comm.close( );
    
    ofstream out_source( source_file.c_str( ) );
    out_source << RUN + "/" + SUBDIR << endl;
    out_source.close( );
  }
  
  // Common part.
  String common = "PRE=" + PRE + " DATA=" + DATA + " RUN=" + RUN;
  String the_command;
  
  // Generate assembly markup info.
  if ( MARKUP ) {
// MakeDepend: dependency RunMarkup
    the_command
      = "RunMarkup " + common
      + " SUBDIR=" + SUBDIR
      + " TAG_REPEATS=True FULL_REPEATS=False MAKE_READABLE=True";
    
    cout << Date( ) << ": running RunMarkup." << endl;
    RunCommandWithLog( the_command, log_file );

    Cp( sub_dir + "/*.xml.gz*", distrib_dir );
    Cp( sub_dir + "/*.gff.gz*", distrib_dir );

// MakeDepend: dependency QualAnnotator
    the_command
      = "QualAnnotator " + common
      + " SUBDIR=" + SUBDIR
      + " RUN_MARKUP_N50=True";
    
    cout << Date( ) << ": running QualAnnotator." << endl;
    RunCommandWithLog( the_command, log_file );

    Cp( sub_dir + "/*.MarkupN50.txt", distrib_dir );
  }

  // Nothing left to do.
  if ( MARKUP_ONLY ) {
    log << Date( ) << ": ForDistributionConst completed successfully" << endl;
    log.close( );
    return 0;
  }
  
  // Some assembly stats.
// MakeDepend: dependency BasicAssemblyStats
  the_command
    = "BasicAssemblyStats " + common
    + " SUBDIR=" + SUBDIR
    + " ARCHIVE=True";

  cout << Date( ) << ": running BasicAssemblyStats." << endl;
  RunCommandWithLog( the_command, log_file );
  Mv( sub_dir + "/" + bastats_file, distrib_dir + "/" + bastats_file );

// MakeDepend: dependency LibStatsOverview
  the_command
    = "LibStatsOverview " + common
    + " SUBDIR=" + SUBDIR
    + " ARCHIVE=True";

  cout << Date( ) << ": running LibStatsOverview." << endl;
  RunCommandWithLog( the_command, log_file );
  Mv( sub_dir + "/" + libstats_file, distrib_dir + "/" + libstats_file );
  
  // Generate csv one-liner with generic info.
// MakeDepend: dependency BasicAssemblyOneLiner
  the_command
    = "BasicAssemblyOneLiner " + common
    + " SUBDIR=" + SUBDIR
    + " TAG_NAME=" + DATA
    + " OUTFILE=" + oneliner_file;
 
  cout << Date( ) << ": running BasicAssemblyOneLiner." << endl;
  RunCommandWithLog( the_command, log_file );
  Mv( sub_dir + "/" + oneliner_file, distrib_dir + "/" + oneliner_file );
  
  // Make readable output.
// MakeDepend: dependency MakeReadableOutput
  the_command
    = "MakeReadableOutput " + common
    + " SUBDIR=" + SUBDIR
    + " OUTDIR=" + BASE_OUT;

  cout << Date( ) << ": running MakeReadableOutput." << endl;
  RunCommandWithLog( the_command, log_file );
  
  // Physical coverage info.
// MakeDepend: dependency PhysicalCoverageByLib
  the_command
    = "PhysicalCoverageByLib " + common
    + " OUTFILE=" + distrib_dir + "/" + physcov_file;

  cout << Date( ) << ": running PhysicalCoverageByLib." << endl;
  RunCommandWithLog( the_command, log_file );
  
  // Generate AGP files.
// MakeDepend: dependency SupersToAgp
  the_command
    = "SupersToAgp " + common
    + " SUBDIR=" + SUBDIR
    + " OUTDIR=" + BASE_OUT;
  
  cout << Date( ) << ": generate AGP file." << endl;
  RunCommandWithLog( the_command, log_file );

  if ( SUPER_TO_FASTA ) {
// MakeDepend: dependency SuperToFasta
    the_command
      = "SuperToFasta " + common
      + " SUBDIR=" + SUBDIR + " QUALS=True";
    
    cout << Date( ) << ": generate super fasta." << endl;
    RunCommandWithLog( the_command, log_file );
  }
  
  // Generate chromosome fasta and qual.
  if ( CHROMOSOME_DATA ) {
    String chr_fasta = "assembly_supers.fasta.gz";
    String chr_qual = "assembly_supers.qual.gz";
    String sub_cfasta_file = sub_dir + "/" + chr_fasta;
    String sub_cqual_file = sub_dir + "/" + chr_qual;
    String cfasta_file = distrib_dir + "/" + chr_fasta;
    String cqual_file = distrib_dir + "/" + chr_qual;
    
    if ( IsRegularFile( sub_agp_file ) )
      Remove( sub_agp_file );
    Cp( agp_file, sub_agp_file );
    
    the_command
// MakeDepend: dependency CreateChromosomeFasta
      = "CreateChromosomeFasta " + common
      + " SUBDIR=" + SUBDIR
      + " OUTDIR=" + SUBDIR
      + " ASSEMBLY_NAME=" + BASE_OUT
      + " AGP_FILE=" + AGPFILE
      + " BASE_OUT=assembly_supers";
    
    RunCommandWithLog( the_command, log_file );
    
// MakeDepend: dependency CreateChromosomeQual
    the_command
      = "CreateChromosomeQual " + common
      + " SUBDIR=" + SUBDIR
      + " OUTDIR=" + SUBDIR
      + " ASSEMBLY_NAME=" + BASE_OUT
      + " AGP_FILE=" + AGPFILE
      + " BASE_OUT=assembly_supers";
    
    RunCommandWithLog( the_command, log_file );

    Mv( sub_cfasta_file, cfasta_file );
    Mv( sub_cqual_file, cqual_file );
  }
  
  // Generate fasta and qual files for unplaced reads.
// MakeDepend: dependency GrabUnplacedReads
  the_command = "GrabUnplacedReads " + common + " SUBDIR=" + SUBDIR;
  
  cout << Date( ) << ": generate unplaced read sequence files." << endl;
  RunCommandWithLog( the_command, log_file );

  String fasta_from = sub_dir + "/unplaced.fasta.gz";
  String fasta_to = distrib_dir + "/unplaced.fasta.gz";
  Cp( fasta_from, fasta_to );

  String qual_from = sub_dir + "/unplaced.qual.gz";
  String qual_to = distrib_dir + "/unplaced.qual.gz";
  Cp( qual_from, qual_to );

  // Copy contigs.fasta in sub_dir.
  String finished_from = data_dir + "/" + KNOWN_CONTIGS;
  String finished_to = sub_dir + "/" + KNOWN_CONTIGS;
  if ( IsRegularFile( finished_from ) )
    Cp( finished_from, finished_to );
  
  // Run EvaluateConsensus and generate Report.
  String workless_common = common;
  if ( workless_common.Contains( "/work" ) )
    workless_common = workless_common.Before( "/work" );
  
// MakeDepend: dependency Assemblez
  the_command
    = "Assemblez " + workless_common
    + " eval_dir=" + SUBDIR
    + " known_contigs=" + KNOWN_CONTIGS
    + " START=EvaluateConsensus";
  
  cout << Date( ) << ": generating Report." << endl;
  RunCommandWithLog( the_command, log_file );

  String report_from = sub_dir + "/" + REPORT_DIR + "/" + REPORT_PRETTY_NAME;
  String report_to = distrib_dir + "/" + REPORT_PRETTY_NAME;
  Cp ( report_from, report_to );

  // Run CalculateReadStats.
  if ( READ_TABLES ) {
// MakeDepend: dependency CalculateReadStats
    the_command
      = "CalculateReadStats " + common
      + " SUBDIR=" + SUBDIR
      + " OUTDIR=" + SUBDIR;
    
    cout << Date() << ": generating Read Statistics Table." << endl;
    RunCommandWithLog( the_command, log_file );
    
    String table_from = sub_dir + "/ReadUsage.Table";
    String table_to = distrib_dir + "/ReadUsage.Table";
    Cp( table_from, table_to );
    
    if ( IsRegularFile(sub_dir + "/Alternate.ReadUsageTable") )
    {
      String alt_from = sub_dir + "/Alternate.ReadUsageTable";
      String alt_to = distrib_dir + "/Alternate.ReadUsageTable";
      Cp( alt_from, alt_to );
    }
  }
  
  // Run CalculateMapStats.
  if ( ASSEMBLY_TABLES ) {
// MakeDepend: dependency CalculateMapStats
    the_command
      = "CalculateMapStats " + common
      + " SUBDIR=" + SUBDIR
      + " OUTFILE=" + sub_dir + "/AssemblyStats.Table"
      + " agpFile=" + agp_file;
    
    cout << Date() << ": generating Assembly Statistics Table." << endl;
    RunCommandWithLog( the_command, log_file );

    String stats_from = sub_dir + "/AssemblyStats.Table";
    String stats_to = distrib_dir + "/AssemblyStats.Table";
    Cp( stats_from, stats_to );
  }

  // Done.
  log << Date( ) << ": ForDistribution completed successfully" << endl;
  log.close( );

}
