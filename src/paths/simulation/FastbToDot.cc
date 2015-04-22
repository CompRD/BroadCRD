///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Basevector.h"
#include "paths/CommonPatherCore.h"
#include "util/RunCommand.h"
// MakeDepend: dependency CommonPather
// MakeDepend: dependency MakeRcDb
// MakeDepend: dependency Unipather
// MakeDepend: dependency UnipathDot

/**
 * FastbToDot
 *
 * Generate the dot file for the HyperKmerPath built from a selected
 * interval from a fastb (for the given K). This is a simple wrapper
 * around UnipathDot.
 *
 * K: kmer size
 * OUTDIR: where all output is saved
 * FASTB: input fastb
 * ID: id in the fastb
 * BEGIN: begin interval
 * LENGTH: length of interval
 */
int main( int argc, char *argv[] )
{
  BeginCommandArguments;
  CommandArgument_Int( K );
  CommandArgument_String( OUTDIR );
  CommandArgument_String( FASTB );
  CommandArgument_Int( ID );
  CommandArgument_Int( BEGIN );
  CommandArgument_Int( LENGTH );
  EndCommandArguments;

  // Dir and file names.
  String ref_head = "target";
  String log_file = OUTDIR + "/FastbToDoc.log";
  String ref_fastb = OUTDIR + "/" + ref_head + ".fastb";

  Mkpath( OUTDIR );
  
  ofstream out( log_file.c_str( ) );
  PrintCommandPretty( out );
  out.close( );

  // Generate target fastb.
  cout << Date( ) << ": generating target fastb" << endl;
  vec<int> select( 1, ID );
  vecbvec all_bases;
  all_bases.SparseRead( FASTB, select, 0 );
  
  vecbvec bases;
  bases.reserve( 1 );
  bases.push_back( basevector( all_bases[ID], BEGIN, LENGTH ) );
  bases.WriteAll( ref_fastb );
  
  // Run CommonPather.
  cout << Date( ) << ": running CommonPather" << endl;
  String theCommand
    = "CommonPather K=" + ToString( K )
    + " READS_IN=\"{" + OUTDIR + "/" + ref_head + ".fastb}\""
    + " PATHS_OUT=\"{" + OUTDIR + "/" + ref_head + ".paths}\"";
  RunCommand( theCommand );

  // Run MakeRcDb.
  theCommand
    = "MakeRcDb K=" + ToString( K )
    + " PRE=" + OUTDIR
    + " DATA= RUN= READS=" + ref_head;
  RunCommand( theCommand );

  // Run Unipather.
  theCommand
    = "Unipather K=" + ToString( K )
    + " PRE=" + OUTDIR
    + " DATA= RUN= READS=" + ref_head;
  RunCommand( theCommand );

  // Run UnipathDot.
  theCommand
    = "UnipathDot K=" + ToString( K )
    + " HEAD=" + OUTDIR + "/" + ref_head;
  RunCommand( theCommand );
  
  // Done.
  cout << Date( ) << ": FastbToDot done" << endl;
  
}
