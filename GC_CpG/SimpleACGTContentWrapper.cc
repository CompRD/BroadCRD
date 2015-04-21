/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Basevector.h"
#include "ReadLocation.h"
#include "util/RunCommand.h"

/**
 * SimpleACGTContentWrapper
 *
 * Run SimpleACGTContent on the following datasets:
 *  1. contigs
 *  2. all reads
 *  3. assembled reads
 *  4. unassembled reads
 *
 * OUTDIR: relative to SUBDIR
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String( SUBDIR );
  CommandArgument_String_OrDefault( OUTDIR, "ACGT_content" );
  EndCommandArguments;

  // Dir and file names.
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/" + SUBDIR;
  String out_dir = sub_dir + "/" + OUTDIR;

  String read_bases = run_dir + "/reads.fastb";

  String contig_bases = sub_dir + "/mergedcontigs.fastb";
  String locs_file = sub_dir + "/mergedcontigs_orig.locs";

  String placed_file = out_dir + "/assembled.ids";
  String unplaced_file = out_dir + "/unassembled.ids";

  Mkdir777( out_dir );

  // Generate assembled/unassembled files.
  cout << Date( ) << ": generating assembled/unassembled lists" << endl;
  vec<Bool> placed( MastervecFileObjectCount( read_bases ), False );
  READ( locs_file, vec<read_location>, locs );
  for (int ii=0; ii<(int)locs.size( ); ii++)
    placed[ locs[ii].ReadId( ) ] = True;
  vec<int> assembled;
  vec<int> unassembled;
  assembled.reserve( placed.size( ) );
  unassembled.reserve( placed.size( ) );
  for (int ii=0; ii<(int)placed.size( ); ii++) {
    if ( placed[ii] )
      assembled.push_back( ii );
    else
      unassembled.push_back( ii );
  }
  WRITE( out_dir + "/assembled.ids", assembled );
  WRITE( out_dir + "/unassembled.ids", unassembled );
  
  // On contigs.
  String theCommand
    = "SimpleACGTContent BASES=" + contig_bases
    + " OUTFILE=" + out_dir + "/contigs.out";
  RunCommand( theCommand );

  // On all reads.
  theCommand
    = "SimpleACGTContent BASES=" + read_bases
    + " OUTFILE=" + out_dir + "/all_reads.out";
  RunCommand( theCommand );
  
  // On assembled reads.
  theCommand
    = "SimpleACGTContent BASES=" + read_bases
    + " SELECT=" + out_dir + "/assembled.ids"
    + " OUTFILE=" + out_dir + "/assembled_reads.out";
  RunCommand( theCommand );
  
  // On unassembled reads.
  theCommand
    = "SimpleACGTContent BASES=" + read_bases
    + " SELECT=" + out_dir + "/unassembled.ids"
    + " OUTFILE=" + out_dir + "/unassembled_reads.out";
  RunCommand( theCommand );

  // Done.
  cout << "\n" << Date( ) << ": all done" << endl;
}
