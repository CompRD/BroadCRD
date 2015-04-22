///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "system/System.h"
#include "util/RunCommand.h"

/**
 * PicardToReads
 *
 * Pre-processing tool to generate the files needed by AssembleViral,
 * from a single picard project.
 *
 * INPUT 1 (it overrides INPUT 2)
 *   <SINGLE_FILE>.{bam,params.txt}
 *
 * INPUT 2
 *   /seq/picard/<FLOWCELL>/picard_subdir/<LANE>/<LIB>/{params.txt,*.bam}
 *   (if PICARD_SUBDIR is given, then picard_subdir is set to it,
 *   otherwise the code expect to find exactly one subdir in
 *   /seq/picard/<FLOWCELL>).
 *
 * OUTPUT:
 *   <OUT_DIR>/<OUT_HEAD>.{fastb,qualb,pairs,mapq,qltout}
 *
 * NB: if the params.txt file does not contain an estimate for the
 *   standard deviation of the lengths of the inserts, this is
 *   approximated as 20% of the insert size.
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;

  CommandArgument_String_OrDefault( SINGLE_FILE, "" );
  
  CommandArgument_String_OrDefault( FLOWCELL, "" );
  CommandArgument_String_OrDefault( LANE, "" );
  CommandArgument_String_OrDefault( LIB, "" );
  CommandArgument_String_OrDefault( PICARD_SUBDIR, "" );

  CommandArgument_String( OUT_DIR );
  CommandArgument_String_OrDefault( OUT_HEAD, "reads" );
  
  EndCommandArguments;

  // Check binaries are in the path.
  if ( ! IsCommandInPath( "samtools" ) ) {
    cout << "Fatal error: samtools is not in your path. Exit.\n" << endl;
    return 1;
  }
  if ( ! IsCommandInPath( "SAM2CRDDump" ) ) {
    cout << "Fatal error: SAM2CRDDump is not in your path. Exit.\n" << endl;
    return 1;
  }

  // Check input files.
  String bam_file;
  String params_file;

  if ( SINGLE_FILE != "" ) {
    bam_file = SINGLE_FILE + ".bam";
    params_file = SINGLE_FILE + ".params.txt";
  }
  else {
    if ( FLOWCELL == "" || LANE == "" || LIB == "" ) {
      if ( ! IsCommandInPath( "samtools" ) ) {
	cout << "Fatal error: if SINGLE_FILE is not defined, then you must\n"
	     << "define FLOWCELL, LANE, and LIB. Exit.\n" << endl;
	return 1;
      }
    }
    
    String in_dir = "/seq/picard/" + FLOWCELL;
    vector<String> allfiles = AllFiles( in_dir );
    vec<String> flow_subdirs;
    if ( PICARD_SUBDIR == "" ) {
      for (size_t ii=0; ii<allfiles.size( ); ii++)
	if ( IsDirectory( in_dir + "/" + allfiles[ii] ) )
	  flow_subdirs.push_back( allfiles[ii] );
      if ( flow_subdirs.size( ) != 1 ) {
	cout << "Fatal error: " << in_dir
	     << " should contain exactly one subdir,\nbut I found "
	     << flow_subdirs.size( ) << " subdirs. Exit.\n" << endl;
	return 1;
      }
    }
    else
      flow_subdirs.push_back( PICARD_SUBDIR );
    in_dir += "/" + flow_subdirs[0] + "/" + LANE + "/" + LIB;
    allfiles = AllFiles( in_dir );
    vec<String> bam_files;
    for (size_t ii=0; ii<allfiles.size( ); ii++) {
      String file_name = in_dir + "/" + allfiles[ii];
      if ( file_name.Contains( ".bam", -1 ) )
	bam_files.push_back( file_name );
    }
    if ( bam_files.size( ) != 1 ) {
      cout << "Fatal error: " << in_dir << " should contain exactly one bam,\n"
	   << "while I found " << bam_files.size( ) << ". Exit.\n" << endl;
      return 1;
    }
    bam_file = bam_files[0];
    params_file = in_dir + "/params.txt";
  }
  
  if ( ! IsRegularFile( params_file ) ) {
    cout << "Fatal error: " << params_file << " not found. Exit.\n" << endl;
    return 1;
  }
  
  // Dir and file names.
  String picard_outdir = OUT_DIR + "/picard";
  String bam_link = picard_outdir + "/orig.bam";
  String params_link = picard_outdir + "/orig.params";

  String out_head = OUT_DIR + "/" + OUT_HEAD;

  // Symlink orig files in pcard_outdir.
  if ( FLOWCELL != "" ) {
    Mkpath( picard_outdir );
    SymlinkForce( bam_file, bam_link );
    SymlinkForce( params_file, params_link );
  }

  // Load insert sep from params.txt
  int insert_size = -666;
  int insert_dev = -666;
  {
    ifstream in( params_file.c_str( ) );
    String line;
    while( in ) {
      getline( in, line );
      if ( line.Contains( "EXPECTED_INSERT_SIZE=" ) ) {
	String size_str = line.After( "EXPECTED_INSERT_SIZE=" );
	if ( size_str.IsInt( ) ) insert_size = size_str.Int( );
      }
      if ( line.Contains( "EXPECTED_INSERT_STDEV=" ) ) {
	String dev_str = line.After( "EXPECTED_INSERT_STDEV=" );
	if ( dev_str.IsInt( ) ) insert_dev = dev_str.Int( );
      }
    }
  }
  if ( insert_size == -666 ) {
    cout << "Warning: no expected insert size found in the params file,\n"
	 << "using 650 as a default.\n" << endl;
    insert_size = 650;
  }
  if ( insert_dev = -666 ) {
    insert_dev = Max( 1, int( 0.2 * double( insert_size ) ) );
  }

  // Pipe SAM2CRDDump to samtools view.
  String theCommand
    = "samtools view -h " + bam_file
    + " | SAM2CRDDump LIBINFO= OUT_HEAD=" + out_head
    + " SEP=" + ToString( insert_size )
    + " DEV=" + ToString( insert_dev );
  RunCommand( theCommand );

  // Done.
  cout << Date( ) << ": PicardToReads done" << endl;
  
}
