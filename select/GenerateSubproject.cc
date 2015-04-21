/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "ReadPairing.h"
#include "TokenizeString.h"
#include "random/Shuffle.h"
#include "system/System.h"
#include "feudal/BinaryStream.h"

/**
 * GenerateSubproject
 *
 * Create a project using a random subset of reads or read pairs (up
 * to the desired coverage) from the project in DATA_IN . Save the
 * output in DATA_OUT.
 *
 * If NAMES_ONLY is set to true, then it saves only the names of the
 * selected reads (use NamesToProject to parse the original input
 * files and generate subproject data).
 *
 * Warning: if reads are not randomized by pairs, it will assume all
 * the reads are unpaired.
 *
 * PAIRS: randomize reads in pairs (if false randomize reads)
 * NAMES_ONLY: save names of selected ids only (see comment above)
 * COVERAGE: desired coverage
 * SEED: for randomization
 * QSEL: use basese of qual >= QSEL to compute coverage
 */
int main( int argc, char *argv[] )
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA_IN );
  CommandArgument_String( DATA_OUT );
  CommandArgument_String_OrDefault( RUN_IN, "run/work" );
  CommandArgument_String_OrDefault( RUN_OUT, "run/work" );
  CommandArgument_Bool_OrDefault( PAIRS, True );
  CommandArgument_Bool_OrDefault( NAMES_ONLY, True );
  CommandArgument_Double_OrDefault( COVERAGE, 10.0 );
  CommandArgument_Int_OrDefault( SEED, 666 );
  CommandArgument_Int_OrDefault( QSEL, 20 );
  EndCommandArguments;
  
  // Input assembly.
  String data_dir_in = PRE + "/" + DATA_IN;
  String run_dir_in = data_dir_in + "/" + RUN_IN;

  String gsize_in = data_dir_in + "/genome.size";
  String config_in = data_dir_in + "/reads_config.xml";
  String organism_in = data_dir_in + "/organism";

  String ids_file_in = run_dir_in + "/reads.ids";
  String pairs_file_in = run_dir_in + "/reads.pairto";
  String bases_file_in = run_dir_in + "/reads.fastb";
  String quals_file_in = run_dir_in + "/reads.qualb";
  String lengths_file_in = run_dir_in + "/reads.lengths";

  // Output assembly.
  String data_dir_out = PRE + "/" + DATA_OUT;
  String run_dir_out = data_dir_out + "/" + RUN_OUT;

  String log_file = data_dir_out + "/GenerateSubproject.log";
  String names_file = data_dir_out + "/select.ids";

  String ids_file_out = run_dir_out + "/reads.ids";
  String bases_file_out = run_dir_out + "/reads.fastb";
  String quals_file_out = run_dir_out + "/reads.qualb";
  String lengths_file_out = run_dir_out + "/reads.lengths";
  
  // Mkdir.
  vec<String> tok;
  vec<char> sep( 1, '/' );
  Tokenize( RUN_OUT, sep, tok );
  String currdir = data_dir_out;
  Mkdir777( currdir );
  for (int ii=0; ii<(int)tok.size( ); ii++) {
    currdir += "/" + tok[ii];
    Mkdir777( currdir );
  }
  
  ofstream log( log_file.c_str( ) );
  PrintCommandPretty( log );
  
  // First copy data_dir files (genome.size, config file, and organism).
  log << Date( ) << ": copying genome size file";
  if ( IsRegularFile( gsize_in ) ) {
    Cp2( gsize_in, data_dir_out );
    log << "\n";
  }
  else
    log << ". Warning! file not found\n";
  
  log << Date( ) << ": copying organism size file";
  if ( IsRegularFile( organism_in ) ) {
    Cp2( organism_in, data_dir_out );
    log << "\n";
  }
  else
    log << ". Warning: file not found\n";

  log << Date( ) << ": copying config file";
  if ( IsRegularFile( config_in ) ) {
    Cp2( config_in, data_dir_out );
    log << "\n";
  }
  else
    log << ". Warning! file not found\n";

  // Load (bases, quals, and pairing info).
  double genome_size = FirstLineOfFile( gsize_in ).Int( );

  log << Date( ) << ": loading bases" << endl;
  vecbvec bases;
  bases.ReadAll( bases_file_in );
  
  log << Date( ) << ": loading quals" << endl;
  vecqvec quals;
  quals.ReadAll( quals_file_in );
  
  log << Date( ) << ": loading pairing info" << endl;
  vec<read_pairing> pairs;
  ReadPairsFile( pairs_file_in, pairs );

  // Evaluate number of bases required.
  longlong target = (longlong)( COVERAGE * genome_size );
  
  // Select reads.
  log << Date ( ) << ": selecting random reads" << endl;

  longlong currcov = 0;
  int npairs = pairs.size( );
  int nreads = bases.size( );

  vec<int> sel_reads;
  vec<int> sel_pairs;
  sel_reads.reserve( bases.size( ) );
  sel_pairs.reserve( pairs.size( ) );
  
  if ( PAIRS ) {
    if ( pairs.size( ) < 1 ) {
      log << "\nFatal error: there are no paired reads.\n" 
	  << "Try running with PAIRS=False.\n"
	  << endl;
      return 0;
    }

    vec<int> randsel;
    Shuffle( npairs, randsel, SEED );

    int endsel = 0;
    while ( currcov < target && endsel < npairs ) {
      const int read_id = randsel[endsel];
      const read_pairing &pair = pairs[read_id];
      sel_pairs.push_back( read_id );
      for (int ii=0; ii<2; ii++) {
	int id = ( ii == 0 ) ? pair.id1 : pair.id2;
	sel_reads.push_back( id );
	const qvec &qual = quals[id];
	for (int jj=0; jj<(int)qual.size( ); jj++)
	  if ( (int)qual[jj] >= QSEL ) currcov++;
      }
      endsel++;
    }
  }
  else {
    vec<int> randsel;
    Shuffle( nreads, randsel, SEED );
    
    int endsel = 0;
    while ( currcov < target && endsel < nreads ) {
      const int read_id = randsel[endsel];
      sel_reads.push_back( read_id );
      const qvec &qual = quals[read_id];
      for (int jj=0; jj<(int)qual.size( ); jj++)
	if ( (int)qual[jj] >= QSEL ) currcov++;
      endsel++;
    }
  }
  
  // Some info.
  double fpairs = 0;
  if ( npairs > 0 )
    fpairs = 100.0 * SafeQuotient( (int)sel_pairs.size( ), npairs );
  double freads = 100.0 * SafeQuotient( (int)sel_reads.size( ), nreads );
  double fcov = SafeQuotient( currcov, (longlong)genome_size );
  String str_fpairs = ToString( fpairs, 2 );
  String str_freads = ToString( freads, 2 );
  String str_cov = ToString( fcov, 2 );
  log << "\n"
      << sel_pairs.size( ) << " pairs were selected, out of an intial set of "
      << npairs << " ("
      << str_fpairs << "%).\n"
      << sel_reads.size( ) << " reads were selected, out of an intial set of "
      << nreads << " ("
      << str_freads << "%).\nThese account for a q"
      << QSEL << " coverage of "
      << str_cov << "x, based on an estimated\ngenome size of "
      << (longlong)genome_size << " (as given in the file genome.size).\n"
      << "\n";
  
  // Save ids.
  log << Date( ) << ": saving read ids" << endl;
  {
    vecString ids_in;
    ids_in.ReadAll( ids_file_in );
    
    // Save selected ids and exit.
    if ( NAMES_ONLY ) {
      ofstream idsout( names_file.c_str( ) );
      for (int ii=0; ii<(int)sel_reads.size( ); ii++)
	idsout << ids_in[ sel_reads[ii] ] << "\n";
      idsout.close( );
      
      log << Date( ) << ": done" << endl;
      log.close( );
      return 0;
    }

    vecString ids_out;
    ids_out.reserve( sel_reads.size( ) );

    for (int ii=0; ii<(int)sel_reads.size( ); ii++)
      ids_out.push_back( ids_in[ sel_reads[ii] ] );
    
    ids_out.WriteAll( ids_file_out );
  }
  
  // Save lengths.
  log << Date( ) << ": saving read lengths" << endl;
  {
    vec<int> lens_out;
    lens_out.reserve( sel_reads.size( ) );
    for (int ii=0; ii<(int)sel_reads.size( ); ii++)
      lens_out.push_back( bases[ sel_reads[ii] ].size( ) );

    BinaryWriter::writeFile( lengths_file_out, lens_out );
  }

  // Save pairs.
  log  << Date( ) << ": saving pairs" << endl;
  {
    vec<read_pairing> pairs_out;
    if ( PAIRS ) {
      pairs_out.reserve( sel_pairs.size( ) );
      for (int ii=0; ii<(int)sel_pairs.size( ); ii++) {
	read_pairing rpair = pairs[ sel_pairs[ii] ];
	rpair.id1 = 2 * ii;
	rpair.id2 = 1 + ( 2 * ii );
	pairs_out.push_back( rpair );
      }
    }
    WritePairs( run_dir_out, pairs_out, sel_reads.size( ) );
  }

  // Save bases.
  log << Date( ) << ": saving bases" << endl;

  int totsel = sel_reads.size( );
  longlong totbases = 0;
  for (int ii=0; ii<(int)sel_reads.size( ); ii++)
    totbases += bases[ sel_reads[ii] ].size( );

  {
    vecbvec bases_out;
    bases_out.Reserve( totbases / 16 + totsel, totsel );
    for (int ii=0; ii<(int)totsel; ii++)
      bases_out.push_back( bases[ sel_reads[ii] ] );
    bases_out.WriteAll( bases_file_out );
  }

  // Save quals.
  log << Date( ) << ": saving quals" << endl;
  {
    vecqvec quals_out;
    quals_out.Reserve( totbases + totsel, totsel );
    for (int ii=0; ii<(int)totsel; ii++)
      quals_out.push_back( quals[ sel_reads[ii] ] );
    quals_out.WriteAll( quals_file_out );
  }
  
  // Done.
  log << Date( ) << ": done" << endl;
  log.close( );
  
}
