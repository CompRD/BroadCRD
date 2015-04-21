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
#include "Bitvector.h"
#include "FeudalMimic.h"
#include "Qualvector.h"
#include "feudal/BinaryStream.h"
#include "kmers/KmersTaggerCore.h"
#include "kmers/ReadPather.h"
#include "kmers/ReadPatherDefs.h"
#include "util/RunCommand.h"

/**
 * KmersTagger
 *
 * Discard low frequency kmers from a set of reads. Output is saved as
 * a fastb, one entry per segment (one read may possibly be broken
 * into multiple segments).
 *
 * Notice that kmer size (K) is hard coded.
 *
 * MIN_FREQ: discard kmers with frequency < MIN_FREQ
 * MIN_READLEN: do not save short segments (kmer length)
 * NUM_THREADS: use all if 0
 * FORCE: do not use cached dictionary
 * VERBOSE: toggle off verbose log
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments; 
  CommandArgument_String( HEAD_IN );
  CommandArgument_String( HEAD_OUT );
  CommandArgument_Int_OrDefault( MIN_FREQ, 2 );
  CommandArgument_Int_OrDefault( MIN_READLEN, 12 );
  CommandArgument_UnsignedInt_OrDefault( NUM_THREADS, 0 );
  CommandArgument_Bool_OrDefault( FORCE, True );
  CommandArgument_Bool_OrDefault( VERBOSE, False );
  CommandArgument_Bool_OrDefault( VALIDATE, False );
  EndCommandArguments;  

  // Constants.
  const unsigned K = 24;

  // Thread control.
  NUM_THREADS = configNumThreads( NUM_THREADS );
  
  // File names.
  String strK = ToString( K );
  String bases_file = HEAD_IN + ".fastb";
  String dict_file = HEAD_IN + ".k" + strK + ".ug.dict";

  String log_file = HEAD_OUT + ".KmersTagger.log";
  String out_bases_file = HEAD_OUT + ".fastb";
  
  vec<String> needed;
  needed.push_back( bases_file );
  
  // Log stream.
  ofstream log( log_file.c_str( ) );
  PrintCommandPretty( log );
  cout << "Sending log to " << log_file << "\n" << endl;
  
  // Load.
  log << Date( ) << ": loading reads" << endl;
  VirtualMasterVec<bvec> bases( bases_file );
  
  KmerDict<K> dictW( 0 );
  KmerDict<K> dictR( 0 );
  bool write_dict = ( FORCE || ! IsRegularFile( dict_file ) );
  if ( write_dict ) {
    log << Date( ) << ": generating (and saving) dict" << endl;
    dictW.process( bases, VALIDATE, NUM_THREADS );
    BinaryWriter::writeFile( dict_file, dictW );
  }
  else {
    log << Date( ) << ": reading dict" << endl;
    BinaryReader::readFile( dict_file, &dictR );
  }
  const KmerDict<K> &dict = write_dict ? dictW : dictR;
  
  // Run KmersTaggerCore.
  vecbvec all_segments;
  KmersTaggerCore<K> ( MIN_FREQ, MIN_READLEN, VERBOSE,
		       dict, bases, all_segments, log );
  
  // Save.
  log << Date( ) << ": saving segments" << endl;
  all_segments.WriteAll( out_bases_file );
  
  // Done.
  String str_done = Date( ) + ": KmersTagger done";
  log << str_done << endl;    
  cout << str_done << endl;    
  log.close( );

}
