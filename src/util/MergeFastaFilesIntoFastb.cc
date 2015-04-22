// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
//

// Merge fasta files into a single fastb file.

#include "MainTools.h"
#include "Basevector.h"
#include "Bitvector.h"
#include "FetchReads.h"
#include "FetchReadsAmb.h"
#include "math/Functions.h"
#include "system/ParsedArgs.h"
#include "ParseSet.h"
#include "STLExtensions.h"
#include "String.h"
#include "system/System.h"
#include "Vec.h"
#include "feudal/IncrementalWriter.h"


/*
 * MergeFastaFilesIntoFastb
 *
 * It converts all the fasta files in a given directory into a single
 * fastb object. It will generate three files: fastb (vecbasevector with
 * bases), fastamb (vecbitvector with ambigouous bases), and ids (vecString
 * with the names of the objects in the fasta files).
 *
 * IN_DIR: full path name for input dir (where fasta files are)
 * OUT_DIR: full path name of output dir
 * TAIL: extensions for fasta input files
 * HEAD: head for output files in OUT_DIR
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( IN_DIR );
  CommandArgument_String( OUT_DIR );
  CommandArgument_String_OrDefault( TAIL, ".fasta" );
  CommandArgument_String_OrDefault( HEAD, "contigs" );
  EndCommandArguments;

  // File names
  String bases_file = OUT_DIR + "/" + HEAD + ".fastb";
  String amb_file = OUT_DIR + "/" + HEAD + ".fastamb";
  String names_file = OUT_DIR + "/" + HEAD + ".ids";

  // Select files.
  vec<String> fasta_files;
  {
    vector<String> all_files = AllFiles( IN_DIR );
    for (int ii=0; ii<(int)all_files.size( ); ii++)
      if ( all_files[ii].Contains( TAIL, -1 ) )
	fasta_files.push_back( all_files[ii] );
    sort( fasta_files.begin( ), fasta_files.end( ) );
  }

  // Prepare to run (remove old runs leftovers).
  Mkdir777( OUT_DIR );

  Remove( bases_file );
  Remove( amb_file );

  vecString contigs_names;
  vecbasevector contigs_bases;
  vecbitvector contigs_amb;

  IncrementalWriter<bitvector> ambWriter(amb_file.c_str());
  IncrementalWriter<bvec> bWriter(bases_file.c_str());

  // Parse fasta files.
  cout << Date( ) << ": parsing " << fasta_files.size( ) << " files:" << endl;

  for (int ii=0; ii<(int)fasta_files.size( ); ii++) {
    Dot( cout, ii );

    const String fasta_file = IN_DIR + "/" + fasta_files[ii];

    // Bases.
    contigs_bases.clear( );
    FetchReads( contigs_bases, 0, fasta_file );
    bWriter.add(contigs_bases.begin(),contigs_bases.end());

    // Ambiguous bases.
    contigs_amb.clear( );
    FetchReadsAmb( contigs_amb, fasta_file );
    ambWriter.add(contigs_amb.begin(),contigs_amb.end());

    // Object names.
    String tmp_file = OUT_DIR + "/tmp";
    String the_comm = "grep \">\" " + fasta_file + " > " + tmp_file;
    System( the_comm );
    ifstream in( tmp_file.c_str( ) );
    while ( in ) {
      String aline;
      getline( in, aline );
      if ( !in ) break;
      contigs_names.push_back( aline.After( ">" ) );
    }
    in.close( );
    Remove( tmp_file );
  }
  cout << "\n";

  // Marge feudal objects and save names.
  cout << Date( ) << ": merging feudal objects and saving names" << endl;
  bWriter.close();
  ambWriter.close();
  contigs_names.WriteAll( names_file );

  // Done.
  cout << Date( ) << ": done" << endl;

}
