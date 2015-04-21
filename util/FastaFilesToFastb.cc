// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

// Collates all fasta files in a given dir into a single fastb file, and
// generates corresponding fastamb, ids files.

// Requirement: each fasta file must contain one and only one fasta entry.

#include "Basevector.h"
#include "Bitvector.h"
#include "MainTools.h"
#include "FetchReads.h"
#include "FetchReadsAmb.h"
#include "system/System.h"



/*
 * FastaFilesToFastb
 *
 * IN_DIR: full path name for input dir (where fasta files are)
 * OUT_DIR: full path name for output
 * OUT_BASE: base name for output files in OUT_DIR
 * TAIL: tail for the files to be selected in WORK_DIR
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( IN_DIR );
  CommandArgument_String( OUT_DIR );
  CommandArgument_String( OUT_BASE );
  CommandArgument_String_OrDefault( TAIL, ".fasta" );
  EndCommandArguments;

  Mkdir777( OUT_DIR );

  // File names.
  String bases_out = OUT_DIR + "/" + OUT_BASE + ".fastb";
  String amb_out = OUT_DIR + "/" + OUT_BASE + ".fastamb";
  String ids_out = OUT_DIR + "/" + OUT_BASE + ".ids";

  // Collect fasta file names.
  vec<String> in_fastas;
  {
    vector<String> all_files = AllFiles( IN_DIR );
    in_fastas.reserve( all_files.size( ) );
    for (int ii=0; ii<(int)all_files.size( ); ii++)
      if ( all_files[ii].Contains( TAIL, -1 ) )
	in_fastas.push_back( all_files[ii] );
    sort( in_fastas.begin( ), in_fastas.end( ) );
  }
  
  // Load one fasta at a time.
  vecbasevector bases;
  vecbitvector amb;
  vec<String> ids;
  cout << Date( ) << ": converting " << in_fastas.size( ) << " fastas" << endl;
  for (int ii=0; ii<(int)in_fastas.size( ); ii++) {
    cout << ii << "\t" << in_fastas[ii] << endl;
    String in_file = IN_DIR + "/" + in_fastas[ii];
    vecbasevector local_bases;
    vecbitvector local_amb;
    FetchReads( local_bases, 0, in_file );
    FetchReadsAmb( local_amb, in_file );
    if ( local_bases.size( ) < 1 ) {
      cout << "  Warning: no bases found" << endl;
      continue;
    }
    if ( local_bases.size( ) > 1 ) {
      cout << "  Fatal error: this file contains " << local_bases.size( )
	   << " entries. Abort." << endl;
      return 0;
    }
    bases.push_back( local_bases[0] );
    amb.push_back( local_amb[0] );
    ids.push_back( in_fastas[ii] );
  }
  
  // Saving.
  cout << Date( ) << ": saving" << endl;
  bases.WriteAll( bases_out );
  amb.WriteAll( amb_out );
  WRITE( ids_out, ids );

  // Done.
  cout << Date( ) << ": done" << endl;
  
}
