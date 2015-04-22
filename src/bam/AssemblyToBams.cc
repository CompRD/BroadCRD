///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "util/RunCommand.h"
#include "bam/AssemblyToBamsCore.h"
#include "bam/AssemblyToBamsDefaults.h"
#include <omp.h>
// MakeDepend: library OMP

/**
 * AssemblyToBams
 *
 * Generate bam files of various AllPaths assembly data (as reads,
 * unibases, contigs, etc). The specified WORK_DIR is used both for
 * input (a default file is loaded in here), and output. Notice the
 * assembly arguments mirror those of RunAllPathsLG.
 *
 * A default file "defaults" must exist in OUT_DIR, and the code loads
 * it to decide what to align. See bam/AssemblyToBamsDefaults.h for
 * details.
 *
 * OUT_DIR: output dir (a valid "defaults" file must be here)
 * REF_HEAD: head of reference genome
 * FORCE: do not used cached aligns
 * NUM_THREADS: use all available if 0
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;

  // These mirror some of the args of RunAllPathsLG.
  CommandArgument_String( PRE );
  CommandArgument_String( REFERENCE_NAME );
  CommandArgument_String( DATA_SUBDIR );
  CommandArgument_String( RUN );

  // Args specific to AssemblyToBams.
  CommandArgument_String( OUT_DIR );
  CommandArgument_String( REF_HEAD );
  CommandArgument_Bool_OrDefault( FORCE, False );
  CommandArgument_Int_OrDefault( NUM_THREADS, 0 );
  
  // Needed binaries / jar files.
  CommandArgument_String_OrDefault( PICARD_DIR, "/seq/software/picard/current/bin/" );
  CommandArgument_String_OrDefault( BWA_BIN, "/seq/software/picard/current/3rd_party/bwa/bwa" );

  EndCommandArguments;
  
  // Dir and file names.
  String ref_dir   = PRE + "/" + REFERENCE_NAME;
  String data_dir  = ref_dir  + "/" + DATA_SUBDIR;
  String run_dir   = data_dir + "/" + RUN;

  String defaults_file = OUT_DIR + "/defaults";
  String ref_fasta = REF_HEAD + ".fasta";
  String ref_lookup = REF_HEAD + ".lookup";
  String ref_dict = REF_HEAD + ".dict";
  
  String create_seq_dict = PICARD_DIR + "/CreateSequenceDictionary.jar";
  
  // Make sure Java-1.6 and samtools are in the path.
  if ( ! CheckJavaVersion( OUT_DIR ) ) return 1;
  TestExecutableByWhich( "samtools" );
  
  // Thread control.
  NUM_THREADS = configNumThreads( NUM_THREADS );
  omp_set_num_threads( NUM_THREADS );

  // Needed files.
  vec<String> needed;
  needed.push_back( defaults_file );
  needed.push_back( ref_fasta );
  needed.push_back( create_seq_dict );
  needed.push_back( BWA_BIN );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;
  
  // Load defaults.
  cout << Date( ) << ": loading defaults" << endl;
  AssemblyToBamsDefaults defaults( defaults_file );
  int nlibs = defaults.NLibraries( );
  
  // Prepare reference.
  PrepareReference( ref_fasta, create_seq_dict, BWA_BIN, OUT_DIR, FORCE );
  
  // Align original reads.
  for (int ii=0; ii<nlibs; ii++) {
    String name = defaults.LibName( ii );
    String in_bam_file = defaults.LibBamFile( ii );

    AlignReadsToReference( in_bam_file, OUT_DIR, name,
			   PICARD_DIR, BWA_BIN, FORCE );
  }
  
  // Align contigs.
  for (int a_id=0; a_id<defaults.NAssemblies( ); a_id++) {
    String head_name = defaults.AssemblyTag( a_id );
    String fastb_file = defaults.AssemblyHead( a_id ) + ".contigs.fastb";

    AlignContigsToReference( fastb_file, OUT_DIR, head_name, FORCE );
  }
  
  // Align unibases.
  for (int u_id=0; u_id<defaults.NUnibases( ); u_id++) {
    String head_name = defaults.UnibasesTag( u_id );
    String fastb_file = defaults.UnibasesFile( u_id );
    bool FW_ONLY = True;
    
    AlignContigsToReference( fastb_file, OUT_DIR, head_name, FORCE, FW_ONLY );
  }

  // Done.
  cout << Date( ) << ": AssemblyToBams done" << endl;
  
}
