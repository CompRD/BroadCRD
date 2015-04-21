///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "TokenizeString.h"
#include "VecUtilities.h"
#include "bam/AssemblyToBamsCore.h"
#include "util/RunCommand.h"

/**
 * AlignBamToFasta
 *
 * Align with bwa a bam file against a reference.
 *
 * REF_FASTA: full path name to reference fasta file
 * BAM_FILE: full path name of input bam file
 * OUT_DIR: where all output will be saved
 * OUT_HEAD: name for output subdir and output bam (defaults to BAM_FILE)
 * PICARD_DIR: directory with picard jar files
 * BWA_DIR: full path name to bwa binary
 * FORCE: force-regenerate all output files
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( REF_FASTA );
  CommandArgument_String( BAM_FILE );
  CommandArgument_String( OUT_DIR );
  CommandArgument_String_OrDefault( OUT_HEAD, "" );
  CommandArgument_String_OrDefault( PICARD_DIR, "/seq/software/picard/current/bin/" );
  CommandArgument_String_OrDefault( BWA_BIN, "/seq/software/picard/current/3rd_party/bwa/bwa" );
  CommandArgument_Bool_OrDefault( FORCE, False );
  EndCommandArguments;

  Mkpath( OUT_DIR );

  // File names.
  String create_seq_dict = PICARD_DIR + "/CreateSequenceDictionary.jar";

  // Check files exist.
  vec<String> needed;
  needed.push_back( REF_FASTA );
  needed.push_back( BAM_FILE );
  needed.push_back( create_seq_dict );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;

  // Make sure Java-1.6 and samtools are in the path.
  if ( ! CheckJavaVersion( OUT_DIR ) ) return 1;
  TestExecutableByWhich( "samtools" );
  
  // Prepare reference.
  PrepareReference( REF_FASTA, create_seq_dict, BWA_BIN, OUT_DIR, FORCE );
  
  // Align reads.
  String name = OUT_HEAD;
  if ( name == "" ) name = Basename( BAM_FILE ).Before( ".bam" );
  AlignReadsToReference( BAM_FILE, OUT_DIR, name, PICARD_DIR, BWA_BIN, FORCE );
  
  // Done
  cout << Date( ) << ": AlignBamToFasta done" << endl;
  
}
