///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "PairsManager.h"
#include "Qualvector.h"
#include "math/Functions.h"
#include "util/RunCommand.h"

/**
 * ReadsToBam
 *
 * Convert fastb/qualb/pairto to bam. This is done in two steps: first
 * convert input data to fastqs, and then fastqs to bam (the latter is
 * done via a system call to a java application).
 *
 * WARNING: the reads are assumed to come from just one library.
 *
 * If LIB_NAME is empty, a default will be generated based on the
 * pairing info, as: lib_SEP_SD, where SEP and SD are the separation
 * and stdev of the library.
 *
 * READS: it will load <READS>.{fastb,qualb,pairto}
 * BASE_OUT: base name for output files (defaulted to READS)
 * LANE_ID: flowcell lane
 * LIB_NAME: name of library (see above)
 * FASTQ_TO_SAM: jar executable name
 * SAMPLE_NAME: optional argument for FASTQ_TO_SAM
 * FASTQ_ONLY: generate fastq files and stop (skip conversion to bam)
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( READS );
  CommandArgument_String_OrDefault( BASE_OUT, "" );
  CommandArgument_String_OrDefault( LANE_ID, "1" );
  CommandArgument_String_OrDefault( LIB_NAME, "" );
  CommandArgument_String_OrDefault( FASTQ_TO_SAM, "/seq/software/picard/current/bin/FastqToSam.jar" );
  CommandArgument_String_OrDefault( SAMPLE_NAME, "na" );
  CommandArgument_Bool_OrDefault( FASTQ_ONLY, False );
  EndCommandArguments;
  
  // Check arguments.
  if ( BASE_OUT == "" ) BASE_OUT = READS;

  // Dir and file names.
  String base_dir = BASE_OUT;
  if ( BASE_OUT.Contains( "/" ) ) base_dir = base_dir.RevBefore( "/" );
  Mkpath( base_dir );
  
  String bases_file = READS + ".fastb";
  String quals_file = READS + ".qualb";
  String pairs_file = READS + ".pairs";

  String log_file = BASE_OUT + ".ReadsToBam.log";
  String fastq1_file = BASE_OUT + "1.fastq";
  String fastq2_file = BASE_OUT + "2.fastq";
  String bam_file = BASE_OUT + ".bam";
  String fastqtobam_log = BASE_OUT + ".FastqToSam.log";

  ofstream log( log_file.c_str( ) );
  PrintCommandPretty( log );

  cout << "Sending log to " << log_file << "\n" << endl;
  
  // Load.
  log << Date( ) << ": loading bases" << endl;
  vecbvec bases( bases_file );

  vecqvec quals;
  if ( IsRegularFile( quals_file ) ) {
    log << Date( ) << ": loading quals" << endl;
    quals.ReadAll( quals_file );
  }
  else {
    log << Date( ) << ": no quals found! Pretending all bases q50" << endl;
  }
  
  log << Date( ) << ": loading pairing info" << endl;
  PairsManager pairs( pairs_file );
  
  // Lib name.
  String lib_name = LIB_NAME;
  if ( lib_name == "" ) {
    if ( pairs.nPairs( ) > 0 ) {
      String str_sep = ToString( pairs.sep( 0 ) );
      String str_sd = ToString( pairs.sd( 0 ) );
      lib_name = "lib_" + str_sep + "_" + str_sd;
    }
  }

  // Generate fastqs.
  log << Date( ) << ": generating fastqs" << endl;
  ofstream out1( fastq1_file.c_str( ) );
  ofstream out2( fastq2_file.c_str( ) );
  
  for (size_t pair_id=0; pair_id<pairs.nPairs( ); pair_id++) {
    longlong id1 = pairs.ID1( pair_id );
    longlong id2 = pairs.ID2( pair_id );

    String str_id1 = ToString( id1 );
    String str_id2 = ToString( id2 );
    String str_lane = ToString( LANE_ID );
    String strId = lib_name + ":" + str_lane + ":" + str_id1 + "_" + str_id2;

    String b1 = bases[id1].ToString( );
    String q1( b1.isize( ) );
    for (size_t jj=0; jj<bases[id1].size( ); jj++) {
      if ( quals.size( ) > 0 ) {
	if ( quals[id1][jj] == 0 ) b1[jj] = 'N';
	q1[jj] = Min( 93, (int)quals[id1][jj] ) + 33;
      }
      else
	q1[jj] = 83;
    }
    out1 << "@" << strId << "/1\n"
	 << b1 << "\n"
	 << "+" << strId << "/1\n"
	 << q1 << "\n";
    
    String b2 = bases[id2].ToString( );
    String q2( b2.isize( ) );
    for (size_t jj=0; jj<bases[id2].size( ); jj++) {
      if ( quals.size( ) > 0 ) {
	if ( quals[id2][jj] == 0 ) b2[jj] = 'N';
	q2[jj] = Min( 93, (int)quals[id2][jj] ) + 33;
      }
      else
	q2[jj] = 83;
    }
    out2 << "@" << strId << "/2\n"
	 << b2 << "\n"
	 << "+" << strId << "/2\n"
	 << q2 << "\n";
  }
  out1.close( );
  out2.close( );

  // Early stop.
  if ( FASTQ_ONLY ) {
    log << Date( ) << ": done" << endl;
    return 0;
  }
  
  // Run FastqToBam.
  log << Date( ) << ": running FastqToSam" << endl;
  String theCommand
    = "java -Xmx8g -jar " + FASTQ_TO_SAM
    + " F1=" + fastq1_file
    + " F2=" + fastq2_file
    + " O=" + bam_file
    + " LB=" + lib_name
    + " TMP_DIR=" + base_dir
    + " SM=" + SAMPLE_NAME
    + " QUALITY_FORMAT=Standard";
  RunCommandWithLog( theCommand, fastqtobam_log );

  // Done
  log << Date( ) << ": done" << endl;
  log.close( );

  // Done.
  cout << Date( ) << ": done" << endl;
  return 0;
  
}
