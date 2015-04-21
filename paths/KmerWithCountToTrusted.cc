/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "kmers/KmerRecord.h"
#include "paths/KmerPath.h"
#include "feudal/BinaryStream.h"

/**
 * FunctionCore
 */
template <int K>
void FunctionCore( const String &prefix,
		   const String &out_tail,
		   const int &min_count );

/**
 * KmersWithCountToTrusted
 *
 * Creates a kmers.trusted file by taking a kmer_with_count input file
 * and accepting only kmers with count >= MIN_COUNT.
 *
 * It needs:
 *   PREFIX (full path name for kmers with count)
 *   PREFIX.paths.kK (paths)
 *   PREFIX.paths_rc.kK (rc paths)
 * 
 * Warning: kmers_with_count and paths must be in sync! You can get
 * this, for example, by running KmerWithCountToVecBVec with
 * MIN_COUNT=0, followed by CommonPather (running with MIN_COUNT > 0
 * would put the files out of sync).
 *
 * K: size of kmers
 * PREFIX: it will look for PREFIX* files (kmer_with_count, paths)
 * MIN_COUNT: min count to tag a kmer as trusted 
 * OUT_TAIL: output will be saved as PREFIX.OUT_TAIL
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_Int( K );
  CommandArgument_String( PREFIX );
  CommandArgument_Int( MIN_COUNT );
  CommandArgument_String_OrDefault( OUT_TAIL, "ids" );
  EndCommandArguments;
  
  // Templatized function call to FunctionCore
#define FUNCTION(K_) FunctionCore<K_> ( PREFIX, OUT_TAIL, MIN_COUNT );  
  DISPATCH_ON_K(K, FUNCTION);
  
  return 0;
}

/**
 * FunctionCore
 */
template <int K>
void FunctionCore( const String &prefix,
		   const String &out_tail,
		   const int &threshold )
{
  cout << Date() << ": Begin." << endl;
  
  // File names.
  String str_k = ".k" + ToString( K );
  String kmers_file = prefix;
  String paths_fw_file = prefix + ".paths" + str_k;
  String paths_rc_file = prefix + ".paths_rc" + str_k;
  String out_file = prefix + "." + out_tail;

  // Check files exist.
  if ( !IsRegularFile( kmers_file ) ) {
    cout << "\nFatal error, " << kmers_file << " not found.\n" << endl;
    return;
  }
  if ( !IsRegularFile( paths_fw_file ) ) {
    cout << "\nFatal error, " << paths_fw_file << " not found.\n" << endl;
    return;
  }
  if ( !IsRegularFile( paths_rc_file ) ) {
    cout << "\nFatal error, " << paths_rc_file << " not found.\n" << endl;
    return;
  }
  
  // Load.
  cout << Date( ) << ": loading kmers_with_count" << endl;
  vec<kmer_with_count<K> > kmers;
  BinaryReader::readFile( kmers_file, &kmers );
  size_t n_kmers = kmers.size();
  
  cout << Date() << ": loading fw paths" << endl;
  vecKmerPath paths( paths_fw_file );

  cout << Date() << ": loading rc paths" << endl;
  vecKmerPath paths_rc( paths_rc_file );
  
  // If these ForceAsserts fail, the kmers_with_count file is not from
  // the same source as the paths_files
  ForceAssert(paths   .size() == n_kmers);
  ForceAssert(paths_rc.size() == n_kmers);
  
  // Analyze all the kmers in the kmers_with_count file.
  cout << "\nParsing " << n_kmers << " kmers (. = 1 million kmers)\n";
  vec<kmer_id_t> trusted_ids;
  for (size_t i = 0; i < n_kmers; i++) {
    if ( i % 1000000 == 0 ) Dot( cout, i / 1000000 );
    int count = kmers[i].Count();

    if (count < threshold) continue;
    
    // Mark this kmer id as trusted.
    trusted_ids.push_back(paths   [i].GetKmer(0));
    trusted_ids.push_back(paths_rc[i].GetKmer(0));
  }
  cout << "\n" << endl;
  
  UniqueSort(trusted_ids);
  
  // Write to output file
  cout << Date( ) << ": saving" << endl;
  BinaryWriter::writeFile( out_file, trusted_ids );
  
  cout << Date( ) << ": done" << endl;
}

/**
 * Instantiate templates
 */
#define INSTANTIATE_FUNCTION_CORE( K_, dummy )   \
template void FunctionCore<K_>( const String &, \
				const String &, \
                                const int &)

FOR_ALL_K( INSTANTIATE_FUNCTION_CORE, unused );
