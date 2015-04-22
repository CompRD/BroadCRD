///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Basevector.h"
#include "Fastavector.h"
#include "util/RunCommand.h"
  
/**
 * ParseIndexedJumpAnnot
 *
 * Parse a .annot file from the IndexedJump analysis, and generate a
 * vectors fasta file (only fw barcodes are saved).
 *
 * INPUT:
 *   <IN_FILE>
 *
 * OUTPUT:
 *   <OUT_DIR>/barcodes{A,B,C,D}.{fasta,fastb,ids}
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( IN_FILE );
  CommandArgument_String( OUT_DIR );
  EndCommandArguments;

  String log_file = OUT_DIR + "/ParseIndexedJumpAnnot.log";

  String barcodesA_ids_file = OUT_DIR + "/barcodesA.ids";
  String barcodesB_ids_file = OUT_DIR + "/barcodesB.ids";
  String barcodesC_ids_file = OUT_DIR + "/barcodesC.ids";
  String barcodesD_ids_file = OUT_DIR + "/barcodesD.ids";
  
  String barcodesA_fasta_file = OUT_DIR + "/barcodesA.fasta";
  String barcodesB_fasta_file = OUT_DIR + "/barcodesB.fasta";
  String barcodesC_fasta_file = OUT_DIR + "/barcodesC.fasta";
  String barcodesD_fasta_file = OUT_DIR + "/barcodesD.fasta";
  
  String barcodesA_fastb_file = OUT_DIR + "/barcodesA.fastb";
  String barcodesB_fastb_file = OUT_DIR + "/barcodesB.fastb";
  String barcodesC_fastb_file = OUT_DIR + "/barcodesC.fastb";
  String barcodesD_fastb_file = OUT_DIR + "/barcodesD.fastb";
  
  vec<String> needed( 1, IN_FILE );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;

  // Log stream.
  ofstream log( log_file.c_str( ) );
  cout << "Sending log to " << log_file << endl;
  PrintCommandPretty( log );
  
  // Streams for fastas.
  ofstream outA( barcodesA_fasta_file.c_str( ) );
  ofstream outB( barcodesB_fasta_file.c_str( ) );
  ofstream outC( barcodesC_fasta_file.c_str( ) );
  ofstream outD( barcodesD_fasta_file.c_str( ) );

  // Parse input file, generate fastas.
  vec<int> n_barcodes( 4, 0 );   // counter for As, Bs, Cs, and Ds
  ifstream in( IN_FILE.c_str( ) );
  while( in ) {
    String line;
    vec<String> tokens;
    getline( in, line );
    if ( !in ) break;
    Tokenize( line, tokens );
    
    if ( tokens.size( ) != 2 ) continue;

    const String &fasta = tokens[0];
    const String &name = tokens[1];
    if ( !( name.Contains( "Idx" ) && name.Contains( "F" ) ) ) continue;

    ofstream *out = 0;
    if ( name.Contains( "FA" ) ) { out = &outA; n_barcodes[0]++; }
    else if ( name.Contains( "FB" ) ) { out = &outB; n_barcodes[1]++; }
    else if ( name.Contains( "FC" ) ) { out = &outC; n_barcodes[2]++; }
    else if ( name.Contains( "FD" ) ) { out = &outD; n_barcodes[3]++; }
    else ForceAssert( 1 == 0 );

    *out << ">" << name.Before( "F" ) << "\n" << fasta << "\n";
  }
  in.close( );

  // Minimal logging.
  log << " A barcodes: " << n_barcodes[0] << "\n"
      << " B barcodes: " << n_barcodes[1] << "\n"
      << " C barcodes: " << n_barcodes[2] << "\n"
      << " D barcodes: " << n_barcodes[3] << "\n"
      << endl;
  
  // Close fasta streams.
  outA.close( );
  outB.close( );
  outC.close( );
  outD.close( );

  // Reload and dump fastb, and ids.
  for (int ii=0; ii<4; ii++) {
    const String *fasta_file = 0;
    const String *fastb_file = 0;
    const String *ids_file = 0;
    if ( ii == 0 ) {
      fasta_file = &barcodesA_fasta_file;
      fastb_file = &barcodesA_fastb_file;
      ids_file = &barcodesA_ids_file;
    }
    else if ( ii == 1 ) {
      fasta_file = &barcodesB_fasta_file;
      fastb_file = &barcodesB_fastb_file;
      ids_file = &barcodesB_ids_file;
    }
    else if ( ii == 2 ) {
      fasta_file = &barcodesC_fasta_file;
      fastb_file = &barcodesC_fastb_file;
      ids_file = &barcodesC_ids_file;
    }
    else if ( ii == 3 ) {
      fasta_file = &barcodesD_fasta_file;
      fastb_file = &barcodesD_fastb_file;
      ids_file = &barcodesD_ids_file;
    }
    else
      ForceAssert( 1 == 0 );

    vec<String> names;
    vec<fastavector> fastas;
    LoadFromFastaFile( *fasta_file, fastas, names );
    
    vecbvec bases;
    bases.reserve( fastas.size( ) );
    for (size_t ii=0; ii<fastas.size( ); ii++)
      bases.push_back( fastas[ii].ToBasevector( ) );
    
    bases.WriteAll( *fastb_file );
    WRITE( *ids_file, names );
  }
  
  // Done.
  log << Date( ) << ": ParseIndexedJumpAnnot done" << endl;
  log.close( );
  
}
