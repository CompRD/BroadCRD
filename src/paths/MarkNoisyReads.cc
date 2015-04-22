///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Basevector.h"
#include "PairsManager.h"
#include "feudal/BinaryStream.h"

/**
 * MarkNoisyReads
 *
 * Find pairs of reads that contain too many A's
 * 
 *   INPUT FILES:
 *      READS_IN.fastb
 *      READS_IN.pairto
 *
 *   OUTPUT FILES:
 *      READS_IN.isnt_noisy
 *  
 */


int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( READS_IN );
  CommandArgument_Double_OrDefault_Doc( MAX_FRAC_A, 0.9, 
      "maximum fraction of A's present in pairs of reads."); 
  CommandArgument_Bool_OrDefault_Doc( REVERSED, False, 
      "reads are revere-complemente of original reads, search for T's."); 
  EndCommandArguments;

  // Check arguments.

  ForceAssert( MAX_FRAC_A <= 1.0 && MAX_FRAC_A >= 0.0 );
  
  // File names.
  String in_reads_file   = READS_IN + ".fastb";
  String in_pairs_file   = READS_IN + ".pairs";
  String out_select_file = READS_IN + ".isnt_noisy";
  
 
  // Reads to be kept.
  longlong nreads = MastervecFileObjectCount( in_reads_file );
  cout << "number of reads in " << in_reads_file << " = " << nreads << endl;
  vec<Bool> reads_to_keep( nreads, True );
  
  // Loading pairs
  cout << Date() << " loading pairs" << endl;
  PairsManager pairs( in_pairs_file );
  size_t npairs = pairs.nPairs();
  size_t nreadsCheck = pairs.nReads();
  cout << "number of reads in " << in_pairs_file << " = " << nreadsCheck << endl;
  
 
  longlong n_noisy_pairs = 0; 
  // Load reads
  cout << Date() << ": loading reads" << endl;
  vecbvec reads( in_reads_file );
  cout << Date() << ": looking  for noisy reads" << endl;

  for ( size_t pi = 0; pi < npairs; pi++ ){
    longlong ri1 = pairs.ID1(pi);
    longlong ri2 = pairs.ID2(pi);
    int Ns = 0;
    if ( ! REVERSED ){
      for ( uint j = 0; j < reads[ri1].size( ); j++ )
	if ( reads[ri1][j] == BASE_A ) ++Ns;
      for ( uint j = 0; j < reads[ri2].size( ); j++ )
	if ( reads[ri2][j] == BASE_A ) ++Ns;
    }else{
      for ( uint j = 0; j < reads[ri1].size( ); j++ )
	if ( reads[ri1][j] == BASE_T ) ++Ns;
      for ( uint j = 0; j < reads[ri2].size( ); j++ )
	if ( reads[ri2][j] == BASE_T ) ++Ns;
    }
    
    if ( double(Ns) / double( reads[ri1].size() + reads[ri2].size() ) >= MAX_FRAC_A ){
      reads_to_keep[ ri1 ] = False;
      reads_to_keep[ ri2 ] = False;
      n_noisy_pairs++; 
    }
    
  } 
   
  cout << "\n"
       << "reads in input:                 " << nreads << "\n"
       << "pairs in input:                 " << npairs << "\n"
       << "noisy pairs:                    " << n_noisy_pairs << "\n"
       << endl;
  

  cout << Date() << " writing out which reads to keep" << endl;
  BinaryWriter::writeFile( out_select_file, reads_to_keep );
  
  // Done.
  cout << Date( ) << ": done" << endl;
  
}
