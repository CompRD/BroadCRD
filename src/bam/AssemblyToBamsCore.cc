///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "String.h"
#include "TokenizeString.h"
#include "bam/AssemblyToBamsCore.h"
#include "lookup/LookAlign.h"
#include "paths/reporting/ReftigUtils.h"
#include "system/TimeUtils.h"
#include "util/RunCommand.h"
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
// MakeDepend: dependency FastbToFinishQualb
// MakeDepend: dependency MakeLookupTable
// MakeDepend: dependency Qltout2BAM

/**
 * CheckPerlVersion
 */
bool CheckPerlVersion( const String &tmp_dir )
{
  if ( System( "perl --version > /dev/null 2>&1" ) != 0 ) {
    cout << "Fatal error: perl not found in path. Exit.\n" << endl;
    return false;
  }

  String perlversion_file = tmp_dir + "/perl_version.out";
  
  String theCommand = "perl --version";
  RunCommandWithLog( theCommand, perlversion_file );
  ifstream in( perlversion_file.c_str( ) );
  String version = "";
  while ( in ) {
    String line;
    getline( in, line );
    if ( ! in ) break;
    if ( ! line.Contains( "This is perl" ) ) continue;
    vec<String> tokens;
    Tokenize( line, tokens );
    for (int ii=0; ii<tokens.isize( ); ii++) {
      if( tokens[ii].Contains( "v", 0 ) ) {
	version = tokens[ii];
	break;
      }
    }
    if ( version != "" ) break;
  }
  if ( IsRegularFile( perlversion_file ) ) Remove( perlversion_file );

  if ( version == "" ) {
    cout << "Fatal error: unknown perl version. Exit.\n" << endl;
    return false;
  }
  vec<String> tokens;
  Tokenize( version, '.', tokens );
  bool ok = ( tokens.size( ) > 1 && tokens[0] == "v5" && tokens[1] == "10" );
  if ( ! ok ) {
    cout << "Fatal error: you need perl v5.10, but I only found " << version
	 << " in your path.\n" << endl;
    return false;
  }

  return true;
}

/**
 * CheckJavaVersion
 */
bool CheckJavaVersion( const String &tmp_dir )
{
  if ( System( "java -version > /dev/null 2>&1" ) != 0 ) {
    cout << "Fatal error: java not found in path. Exit.\n" << endl;
    return false;
  }

  String javaversion_file = tmp_dir + "/java_version.out";
  
  String theCommand = "java -version";
  RunCommandWithLog( theCommand, javaversion_file );
  ifstream in( javaversion_file.c_str( ) );
  String version = "";
  while ( in ) {
    String line;
    getline( in, line );
    if ( ! in ) break;
    if ( ! line.Contains( "java version" ) ) continue;
    vec<String> tokens;
    Tokenize( line, tokens );
    for (int ii=0; ii<tokens.isize( ); ii++) {
      if( tokens[ii].Contains( "\"", 0 ) ) {
	version = tokens[ii];
	break;
      }
    }
    if ( version != "" ) break;
  }
  if ( IsRegularFile( javaversion_file ) ) Remove( javaversion_file );
  
  if ( version == "" || ! version.Contains( "\"" ) ) {
    cout << "Fatal error: unknown java version. Exit.\n" << endl;
    return false;
  }
  version = version.After( "\"" );
  if ( version.Contains( "\"" ) ) version = version.Before( "\"" );
  vec<String> tokens;
  Tokenize( version, '.', tokens );
  bool ok = ( tokens.size( ) > 1 && tokens[0] == "1" && tokens[1] == "6" );
  if ( ! ok ) {
    cout << "Fatal error: you need java 1.6, but I only found" << version
	 << " in your path.\n" << endl;
    return false;
  }

  return true;
}

/**
 * PrepareReference
 */
void PrepareReference( const String &ref_fasta,
		       const String &createSeqDict_jar,
		       const String &bwa_bin,
		       const String &base_dir,
		       const bool FORCE )
{
  // Dir and file names.
  String out_dir = base_dir + "/reference";
  String out_head = out_dir + "/genome";
  
  String out_fasta = out_head + ".fasta";
  String out_bwt = out_head + ".fasta.bwt";
  String out_dict = out_head + ".dict";
  String out_dict_tmp = out_head + ".dict_tmp";
  String out_lookup = out_head + ".lookup";
  
  Mkpath( out_dir );

  // Link fasta in out_dir.
  if ( FORCE || ! IsRegularFile( out_fasta ) )
    SymlinkForce( ref_fasta, out_fasta );
  
  // Create dictionary.
  if ( FORCE || ! IsRegularFile( out_dict ) ) {
    String theCommand
      = "java -Xmx6g -jar " + createSeqDict_jar
      + " R=" + out_fasta
      + " O=" + out_dict_tmp
      + " TMP_DIR=" + out_dir;
    RunCommand( theCommand );
  }

  // This complicate hack is due to the fact that Qltout2SAM does not
  //  know what to do with the "@HD" line. The correct thing to do
  //  would be to replace the call to createSeqDict_jar with Ted's
  //  RefDesc class, which can save a dictionary optionally skipping
  //  the "@HD" line.
  ifstream tmp_in( out_dict_tmp.c_str( ) );
  ofstream tmp_out( out_dict.c_str( ) );
  String aline;
  while ( tmp_in ) {
    getline( tmp_in, aline );
    if ( ! tmp_in ) break;
    if ( aline.Contains( "@HD", 0 ) ) continue;
    tmp_out << aline << "\n";
  }
  tmp_out.close( );
  if ( IsRegularFile( out_dict_tmp ) ) Remove( out_dict_tmp );

  // Generate index.
  if ( FORCE || ! IsRegularFile( out_bwt ) ) {
    String theCommand
      = bwa_bin + " index -a is " + out_fasta;
    RunCommand( theCommand );
  }

  // Generate lookup table.
  if ( FORCE || ! IsRegularFile( out_lookup ) ) {
    String theCommand
      = "MakeLookupTable SOURCE=" + out_fasta
      + " OUT_HEAD=" + out_head;
    RunCommand( theCommand );
  }

}

/**
 * AlignReadsToReference
 *
 * Use bwa to align reads, and generate .bam and .bai files.
 *
 * in_bam_file: input bam file
 * base_dir, head_name: define output dir (<base_dir>/<head_name>)
 * picard_dir: where picard tools are
 * bwa_bin: full path name to bwa binary
 * FORCE: force-regenerate all output files
 */
void AlignReadsToReference( const String &in_bam_file,
			    const String &base_dir,
			    const String &head_name,
			    const String &picard_dir,
			    const String &bwa_bin,
			    const bool FORCE )
{
  // Dir and file names.
  String out_dir = base_dir + "/" + head_name;

  String ref_fasta = base_dir + "/reference/genome.fasta";

  String bam_orig = out_dir + "/orig.bam";
  String bam_base = Basename( in_bam_file );

  String bam_unmapped = out_dir + "/unmapped.bam";
  String sai_aligned1 = out_dir + "/aligned1.sai";
  String sai_aligned2 = out_dir + "/aligned2.sai";
  String sai_log1 = out_dir + "/aligned1.sai.log";
  String sai_log2 = out_dir + "/aligned2.sai.log";
  String sam_aligned = out_dir + "/aligned.sam";
  String bam_unmarked = out_dir + "/unmarked.bam";
  String metrics_unmarked = out_dir + "/unmarked.summary_metrics";
  String bam_out = out_dir + "/" + head_name + ".bam";
  String metrics_out = out_dir + "/" + head_name + ".metrics";
  
  String revert_sam = picard_dir + "/RevertSam.jar";
  String merge_bam_alignment = picard_dir + "/MergeBamAlignment.jar";
  String collect_metrics = picard_dir + "/CollectAlignmentSummaryMetrics.jar";
  String mark_duplicates = picard_dir + "/MarkDuplicates.jar";

  vec<String> needed;
  needed.push_back( in_bam_file );
  needed.push_back( revert_sam );
  needed.push_back( merge_bam_alignment );
  needed.push_back( collect_metrics );
  if ( ! CheckFilesExist( needed, &cout ) ) return;
  
  Mkpath( out_dir );

  // Nothing to do.
  if ( IsRegularFile( bam_out ) && ! FORCE ) return;

  // Link in original bam file.
  if ( FORCE || ! IsRegularFile( bam_orig ) )
    SymlinkForce( in_bam_file, bam_orig );

  // RevertSam. NB: do not use TMP_DIR with RevertSam! There is a bug
  //  that causes the java job to hang without producing neither any
  //  output, nor a traceback (sante - 2012.02.13).
  if ( FORCE || ! IsRegularFile( bam_unmapped ) ) {
    String args1 = "OQ=true SO=queryname";
    String args2 = "REMOVE_DUPLICATE_INFORMATION=true";
    String args3  ="REMOVE_ALIGNMENT_INFORMATION=true";
    String theCommand
      = "java -Xmx6g -jar " + revert_sam
      + " I=" + bam_orig
      + " O=" + bam_unmapped
      + " " + args1
      + " " + args2
      + " " + args3;
    RunCommand( theCommand );
  }

  // Align reads (paired reads are aligned in two separate groups).
  #pragma omp parallel for
  for (int ii=0; ii<2; ii++) {
    String sai_file = ( ii == 0 ) ? sai_aligned1 : sai_aligned2;
    String log_file = ( ii == 0 ) ? sai_log1 : sai_log2;
    String str_b = ( ii == 0 ) ? "-b1" : "-b2";

    #pragma omp critical
    {
      cout << "Aligning " << str_b << " ends.\n"
	   << "See log file: " << log_file << "\n"
	   << endl;
    }
    
    if ( FORCE || ! IsRegularFile( sai_file ) ) {
      String theCommand
	= bwa_bin + " aln " + ref_fasta
	+ " -q 5 -l 32 -k 2 -t 4 -o 1 -f " + sai_file
	+ " " + str_b
	+ " " + bam_unmapped;
      RunCommandWithLog( theCommand, log_file );
    }
  } 

  // Merge reads from the two groups into a single bam file.
  if ( FORCE || ! IsRegularFile( sam_aligned ) ) {
    String theCommand
      = bwa_bin + " sampe -a 100000 -f "
      + sam_aligned + " "
      + ref_fasta + " "
      + sai_aligned1 + " "
      + sai_aligned2 + " "
      + bam_unmapped + " "
      + bam_unmapped;
    RunCommand( theCommand );
  }
  
  // Run MergeBamAlignment.
  if ( FORCE || ! IsRegularFile( bam_unmarked ) ) {
    String args1 = "PROGRAM_RECORD_ID=BWA PROGRAM_GROUP_VERSION=0.5.7";
    String args2 = "PROGRAM_GROUP_COMMAND_LINE=tk VALIDATION_STRINGENCY=SILENT";
    String args3 = "CLIP_ADAPTERS=false PE=true";
    String theCommand
      = "java -Xmx6g -jar " + merge_bam_alignment
      + " R=" + ref_fasta
      + " UNMAPPED_BAM=" + bam_unmapped
      + " ALIGNED_BAM=" + sam_aligned
      + " OUTPUT=" + bam_unmarked
      + " TMP_DIR=" + out_dir
      + " " + args1
      + " " + args2
      + " " + args3;
    RunCommand( theCommand );
  }

  // Generate summary metrics.
  if ( FORCE || ! IsRegularFile( metrics_unmarked ) ) {
    String theCommand
      = "java -Xmx6g -jar " + collect_metrics
      + " R=" + ref_fasta
      + " I=" + bam_unmarked
      + " O=" + metrics_unmarked
      + " TMP_DIR=" + out_dir;
    RunCommand( theCommand );
  }

  // Mark duplicates (final aligns file).
  String theCommand
    = "java -Xmx6g -jar " + mark_duplicates
    + " I=" + bam_unmarked
    + " O=" + bam_out
    + " M=" + metrics_out
    + " TMP_DIR=" + out_dir;
  RunCommand( theCommand );

  // Generate index.
  theCommand
    = "samtools index " + bam_out;
  RunCommand( theCommand );

  // Clean up.
  if ( IsRegularFile( bam_out ) ) {
    if ( IsRegularFile( bam_unmapped ) ) Remove( bam_unmapped );
    if ( IsRegularFile( sai_aligned1 ) ) Remove( sai_aligned1 );
    if ( IsRegularFile( sai_aligned2 ) ) Remove( sai_aligned2 );
    if ( IsRegularFile( sam_aligned ) ) Remove( sam_aligned );
  }

}

/**
 * AlignContigsToReference
 */
void AlignContigsToReference( const String &in_fastb_file,
			      const String &base_dir,
			      const String &head_name,
			      const bool FORCE,
			      const bool FW_ONLY )
{
  // Dir and file names
  String out_dir = base_dir + "/" + head_name;

  String cgs_file = in_fastb_file;
  String ref_head = base_dir + "/reference/genome";
  String ref_lookup = ref_head + ".lookup";

  String out_head             = out_dir + "/" + head_name;
  String workdir_cgs_file     = out_head + ".fastb";
  String workdir_cgquals_file = out_head + ".qualb";
  String alsA_file            = out_head + ".all.qltout";
  String als_file             = out_head + ".qltout";
  String metrics_file         = out_head + ".metrics";
  String bam_file             = out_head + ".bam";
  String index_file           = out_head + ".bam.bai";
  
  if ( ! IsRegularFile( cgs_file ) ) {
    cout << "\nFATAL ERROR: " << cgs_file << " not found.\nExit." << endl;
    ForceAssert( 1 == 0 );
  }
  
  Mkpath( out_dir );
  
  vec<look_align> hits;
  GetAlignsFast( 12, cgs_file, ref_lookup, als_file, hits, ! FORCE, out_dir );

  if ( FW_ONLY && ( FORCE || ! IsRegularFile( alsA_file ) ) ) {
    Mv( als_file, alsA_file );
    ofstream out( als_file.c_str( ) );
    for (size_t ii=0; ii<hits.size( ); ii++)
      if ( hits[ii].IsQueryFW( ) )
	hits[ii].PrintParseable( out );
    out.close( );
  }
  
  if ( FORCE || ! IsRegularFile( workdir_cgquals_file ) ) {
    SymlinkForce( cgs_file, workdir_cgs_file );
    String theCommand
      = "FastbToFinishQualb BASE_DIR=" + out_dir
      + " HEAD=" + head_name;
    RunCommand( theCommand );
  }
  
  if ( FORCE || ! IsRegularFile( metrics_file ) ) {
    DateString date;
    ofstream mout( metrics_file.c_str( ) );
    mout << "DB_EXTERNAL_ID=sample_" + head_name + "\n"
	 << "DB_LIBRARY_NAME=" + head_name + "\n"
	 << "DB_LANE=none.0\n"
	 << "DB_RUN_NAME=" << date << "_" + head_name + "\n"
	 << "alignment_command=qlt\n";
    mout.close( );
  }
  
  if ( FORCE || ! IsRegularFile( bam_file ) ) {
    String theCommand
      = "Qltout2BAM REF_HEAD=" + ref_head
      + " HEAD_IN=" + out_head
      + " FORCE=True";   // notice FORCE is always true here
    RunCommand( theCommand );
  }
  
  if ( FORCE || ! IsRegularFile( index_file ) ) {
    String theCommand
      = "samtools index " + bam_file;
    RunCommand( theCommand );
  }

}

