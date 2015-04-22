///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "Validates input data for ALLPATHS-LG. Checks that supplied libraries fit "
  "our requirements.";
  
#include "MainTools.h"
#include "String.h"
#include "Vec.h"
#include "Basevector.h"
#include "feudal/BinaryStream.h"
#include "feudal/VirtualMasterVec.h"
#include "PairsManager.h"


enum LibType {frag,jump,long_jump,long_read};

typedef VirtualMasterVec<basevector> VVecBVec;
  
void ErrorMessage(const String message) {
  cout << "[ERROR] " << message << endl; 
}

bool IsFile(const String filename) {
  bool found = IsRegularFile(filename);
  if (!found)
    ErrorMessage("Cannot find file: " + filename);
  return found;
}

bool ErrorIfTrue(bool test, const String message) {
  if (test)
    ErrorMessage(message);
  return test;
}

bool ErrorIfFalse(bool test, const String message) {
  return ErrorIfTrue(!test, message);
}


PM_LibraryStats 
GetUnpairedLibraryStats( const String bases_filename )  {
  PM_LibraryStats stats;
  
  // populate basic library stats
  stats.name = "Unpaired";

  // Use VirtualMasterVec to avoid loading in vecbasevector
  VVecBVec bases((bases_filename).c_str());

  // compute read stats for each library
  for (VVecBVec::const_iterator basesItr = bases.begin(); basesItr != bases.end(); basesItr++) {
    uint32_t read_len= basesItr->size( );
    if (stats.n_reads == 0) {
      stats.min_len = stats.max_len = stats.n_bases = read_len;
    } else {
      if (read_len < stats.min_len) stats.min_len = read_len;
      if (read_len > stats.max_len) stats.max_len = read_len;
      stats.n_bases += read_len;
    }
    stats.n_reads++;
  }
 
  // Compute additional library stats
  stats.mean_len = (stats.n_reads == 0 ? 0 : stats.n_bases / stats.n_reads);
  
  return stats;
}

// Performs various tests on a set of libraries and returns true if no problems were found

bool ValidateReference(const String file_head) {
  String filename = file_head + ".fastb";

  if (!IsFile(filename) )
    return false;

  VVecBVec ref(filename.c_str());

  if (ErrorIfTrue(ref.size() == 0, "No contigs found.") )
    return false;
  
  uint32_t ref_bases = 0;
  for (VVecBVec::const_iterator it = ref.begin(); 
       it != ref.end(); it++) {
    if (ErrorIfTrue(it->size() == 0, "Zero length contig found.") )
      return false;
    ref_bases += it->size(); 
  }

  cout << "Contigs     : " << ref.size() << endl;
  cout << "Length      : " << ref_bases << endl;

  return true;
}


// Performs various tests on a set of libraries and returns true if no problems were found

bool ValidateLibraries(const LibType lib_type, const String file_head, const int K, ostream &log ) {
  bool error = false;

  error |= !IsFile(file_head + ".fastb");

  if (lib_type != long_read) {
    error |= !IsFile(file_head + ".qualb");
    error |= !IsFile(file_head + ".pairs"); 
  }

  if (error) // Missing files, unable to continue
    return false;

  size_t npairs, nlibs, nreads;
  vec<PM_LibraryStats> stats;

  if (lib_type != long_read) {   // Paired libraries
    PairsManager pairs;
    pairs.Read(file_head + ".pairs");
  
    npairs = pairs.nPairs();
    nlibs = pairs.nLibraries();
    nreads = pairs.nReads();
  
    if (ErrorIfTrue( nlibs == 0, "No libraries founds.") )  // Empty files, unable to continue
      return false;

    stats = pairs.getLibraryStats(file_head + ".fastb");

  } else {  // Special case for unpaired libraries

    stats.push_back(GetUnpairedLibraryStats(file_head + ".fastb"));
    npairs = 0;
    nlibs = 1;
    nreads = stats[0].n_reads;
  }

  // Display Library Stats

  cout << "Libraries     : " << nlibs << endl;
  cout << "Pairs         : " << npairs << endl;
  cout << "Total reads   : " << nreads << endl;
  cout << "Paired reads  : " << npairs * 2 << endl;
  cout << "Unpaired reads: " << nreads - (npairs * 2) << endl;
  cout << endl;

  writeLibraryStats( cout, stats );

  // Save output as a csv table.

  String str_type;
  if ( lib_type == frag )      str_type = "frag";
  else if ( lib_type == jump ) str_type = "jump";
  else                         str_type = "long_jump";

  log << "lib_type,lib_name,tot_bases,min_read_len,max_read_len,mean_read_len\n";
  for (size_t lib_id=0; lib_id<stats.size( ); lib_id++) {
    const PM_LibraryStats &lib = stats[lib_id];
    log << str_type << ","
	<< lib.name << ","
	<< lib.n_bases << ","
	<< lib.min_len << ","
	<< lib.max_len << ","
	<< lib.mean_len << "\n";
  }
  log << endl;

  // Basic consistency tests
  
  if (lib_type != long_read) {
    error |= ErrorIfTrue(MastervecFileObjectCount(file_head + ".fastb") != nreads,
			 "Inconsistency found between pairs manager and fastb");

    error |= ErrorIfTrue(MastervecFileObjectCount(file_head + ".qualb") != nreads,
			 "Inconsistency found between pairs manager and qualb");
  }

  // General per library tests
  for ( size_t i = 0; i < nlibs; i++ ) {
    String libname =  stats[i].name;
    if (lib_type == frag) {
      error |= ErrorIfTrue(2*stats[i].max_len + stats[i].sep + 2 * stats[i].sd <= static_cast<unsigned>(K),
			   "Library " + libname + " contains no reads that overlap to create a super read larger K (" + ToString(K) + ").") ;
    }
  }
  
  // Fragment library tests
  if (lib_type == frag) {
    bool found_good_sep = false;
    for ( size_t i = 0; i < nlibs; i++ ) 
      found_good_sep |= (stats[i].sep - 2 * stats[i].sd <= 0 );
    error |= ErrorIfFalse(found_good_sep, "Could not find a fragment library whose reads overlapped.");
  }

  // Jumping library tests
  if (lib_type == jump) {
    bool found_good_sep = false;
    for ( size_t i = 0; i < nlibs; i++ ) 
      found_good_sep |= (stats[i].sep > 500 && stats[i].sep < 10000 );
    error |= ErrorIfFalse(found_good_sep, "Could not find a jumping library with a separation > 1000 and < 10000 bases.");
  }

  // Long jumping library tests
  if (lib_type == long_jump) {
    bool found_good_sep = false;
    for ( size_t i = 0; i < nlibs; i++ ) 
      found_good_sep |= (stats[i].sep > 5000 );
    error |= ErrorIfFalse(found_good_sep, "Could not find a long jumping library with a separation > 5000.");
  }

  // Long read tests
  if (lib_type == long_read) {
    error |= ErrorIfTrue(stats[0].mean_len < static_cast<uint32_t>(K), 
			 "Mean read length < " + ToString(K) + " (CLR_KOUT)");
  }

  return !error;
}





int main( int argc, char *argv[] )
{
  RunTime( );
  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_Int(K);
  CommandArgument_Int(K_LONG);
  CommandArgument_Bool_OrDefault_Doc( WARN_ONLY, False,
    "Don't stop executation of pipeline if a problem is found.")
  CommandArgument_Bool_OrDefault_Doc( LONG_JUMPS, False,
    "Validate long jumping read libraries.")
  CommandArgument_Bool_OrDefault_Doc( LONG_READS, False,
    "Validate long jumping read libraries.")
  CommandArgument_Bool_OrDefault_Doc( REFERENCE, False,
    "Validate optional reference genome.")
  CommandArgument_String_OrDefault_Doc( REPORT, "ValidateAllPathsInputs_core_stats.out",
    "Save to file core stats (used by reporting code).")
  EndCommandArguments;

  String data_dir = PRE + "/" + DATA;

  String frag_reads_head = data_dir + "/frag_reads_orig";
  String jump_reads_head = data_dir + "/jump_reads_orig";
  String long_jump_reads_head = data_dir + "/long_jump_reads_orig";
  String long_reads_head = data_dir + "/long_reads_orig";
  String reference_head = data_dir + "/genome";
  String core_stats_log_file = data_dir + "/" + REPORT;

  ofstream log( core_stats_log_file.c_str( ) );
  PrintCommandPretty( log );
  
  bool error = false;

  // Validate fragment reads

  cout << "Validating Fragment Libraries" << endl
       << "=============================" << endl << endl;

  error |= !ValidateLibraries(frag, frag_reads_head, K, log); 


  cout << endl
       << "Validating Jumping Libraries" << endl
       << "============================" << endl << endl;

  error |= !ValidateLibraries(jump, jump_reads_head, K, log); 


  if (LONG_JUMPS) {
    cout << endl
	 << "Validating Long Jumping Libraries" << endl
	 << "=================================" << endl << endl;
    
    error |= !ValidateLibraries(long_jump, long_jump_reads_head, K, log); 
  }


  if (LONG_READS) {
    cout << endl
	 << "Validating Long Reads" << endl
	 << "=====================" << endl << endl;

    error |= !ValidateLibraries(long_read, long_reads_head, K_LONG, log); 
  }

  if (REFERENCE) {
    cout << endl
	 << "Validating Reference Genome" << endl
	 << "===========================" << endl << endl;

    error |= !ValidateReference(reference_head); 
  }

  cout << endl;
  log.close( );
  
  if (!error)
    cout << "Validation succeeded." << endl;
  if (error && WARN_ONLY) 
    cout << "Validation failed, but proceeding anyway." << endl;
  else if (error) {
    cout << "Validation failed. Aborting RunAllPathsLG." << endl
	 << "Please fix problems noted above and try again." << endl
	 << "For library requirements see the ALLPATHS-LG manual." << endl
	 << "To continue despite the validation failure, use VAPI_WARN_ONLY=True" << endl;
    exit(1);
  }
}
