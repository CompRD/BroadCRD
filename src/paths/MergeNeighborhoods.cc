///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h" // HyperBasevector
#include "paths/HyperKmerPath.h" // HyperKmerPath
#include "paths/InternalMerge.h" // GroupedInternalMergeLG
#include "paths/KmerBaseBroker.h" // KmerBaseBroker
#include "paths/PdfEntry.h" // pdf_entry
#include "paths/ReadsToPathsCoreX.h" // ReadsToPathsCoreY
#include "feudal/BinaryStream.h"


/* MergeNeighborhoods
 *
 * Merge the local-neighborhood assemblies that have been created with
 * LocalizeReadsLG.
 *
 * The central algorithm here is GroupedInternalMergeLG, which performs the
 * actual merge.
 *
 * Josh Burton
 * December 2008
 *
 ******************************************************************************/











int main( int argc, char *argv[] )
{
  RunTime( );
  
  
  // Command-line arguments
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String_OrDefault( READS, "reads" );
  CommandArgument_Int( K );
  
  CommandArgument_Int_OrDefault_Doc( NUM_THREADS, 1, "Number of parallel threads for pathing and KmerParcels" );
  
  // Heuristic parameters for InternalMerge
  CommandArgument_Int_OrDefault_Doc( MIN_ALIGN_LENGTH, 500, "Mark HyperKmerPaths as adjacent if they share a MIN_ALIGN_LENGTH-mer" );
  CommandArgument_Int_OrDefault_Doc( MAX_KMER_FREQ, 7, "Ignore MIN_ALIGN_LENGTH-mers with frequency above this" );
  CommandArgument_Int_OrDefault( MAX_GROUP_SIZE, 3 );
  CommandArgument_Int_OrDefault( MIN_OVERLAP, 10000 );
  CommandArgument_Int_OrDefault( MIN_PROPER_OVERLAP, 10000 );
  CommandArgument_Int_OrDefault( MIN_PROPER_OVERLAP_FINAL, 5000 );
  CommandArgument_Int_OrDefault( MAX_EXTENSION_QUEUE_SIZE, 20 );

  CommandArgument_LongLong_OrDefault( CHECKPOINT_INTERVAL, 4L*60L*60L );
  CommandArgument_Bool_OrDefault( SKIP_FINAL_MERGE, false );
  CommandArgument_Bool_OrDefault( USE_OLD_SEEDS, false );
  EndCommandArguments;
  
  
  
  
  cout << Date( ) << ": Beginning MergeNeighborhoods..." << endl;
  
  
  /*****************************************************************************
   *
   *        LOAD INPUT FILES
   *
   ****************************************************************************/
  
  cout << Date( ) << ": Loading input files" << endl;
  
  // Filenames
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String file_head = run_dir + "/" + READS;
  String kK = ".k" + ToString( K );
  String data_dir = PRE + "/" + DATA;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  
  // Ploidy
  int ploidy = FirstLineOfFile( data_dir + "/ploidy" ).Int( );
  
  // Unipaths
  vecbasevector unibases( file_head + ".unibases" + kK );
  
  // Unipath copy numbers (predicted)
  VecPdfEntryVec CNs( (file_head + ".unipaths.predicted_count" + kK).c_str() );
  
  // HyperBasevectors representing local assemblies
  cout << Date( ) << ": Loading HyperBasevectors from files in "
       << sub_dir << "..." << endl;
  
  int n_seeds = FirstLineOfFile( sub_dir + "/seeds.ids" ).Int( );
  vec<HyperBasevector> local_HBVs( n_seeds );
  for ( int s = 0; s < n_seeds; s++ ) {
    // If any HyperBasevector files are missing, print a warning.
    String HBV_file;
    if (USE_OLD_SEEDS)
      HBV_file = sub_dir + "/seed/" + GetNestedDirsFromKey( s ) + "/hbv";
    else
      HBV_file = sub_dir + "/seed/" + GetDirFromKey( s, 1000 ) + "/" + ToString(s) + ".hbv";
    if ( IsRegularFile( HBV_file ) )
      BinaryReader::readFile( HBV_file, &local_HBVs[s] );
    else
      cout << "WARNING: Could not find HyperBasevector for seed #" << s << endl;
  }
  
  
  /*****************************************************************************
   *
   *        DATASET PRE-PROCESSING: Build some useful auxiliary data structures
   *
   ****************************************************************************/
  
  cout << Date( ) << ": Pre-processing dataset" << endl;
  
  const size_t n_unipaths = unibases.size( );
  
  
  // Read the CNs and find the predicted copy number (CN) of each unipath
  vec<int> predicted_CNs( n_unipaths, -1 );
  for ( size_t i = 0; i < n_unipaths; i++ )
    GetMostLikelyValue( predicted_CNs[i], CNs[i] );
  
  
  
  /*****************************************************************************
   *
   *        MERGE LOCAL ASSEMBLIES INTO GLOBAL ASSEMBLY
   *  This code is adapted from the module LocalizeReadsTail, which is part of
   *  the original RunAllPaths pipeline.
   *
   ****************************************************************************/
  
  
  // Make a vecbasevector out of all of the basevectors in all of the local
  // HyperBasevector objects.
  cout << Date( ) << ": Creating a vecbasevector for the entire assembly" << endl;
  vecbasevector HKP_bases;
  longlong nbases = 0, n_HKP_paths = 0;
  vec<int> base_to_HKP_ID; // maps IDs in 'bases' to HyperKmerPath IDs
  // Pass 1: reserve memory; Pass 2: create the vecbasevector
  for ( int pass = 1; pass <= 2; pass++ ) {
    
    for ( int i = 0; i < n_seeds; i++ ) {
      if ( pass == 1 ) n_HKP_paths += local_HBVs[i].EdgeObjectCount( );
      
      for ( int j = 0; j < local_HBVs[i].EdgeObjectCount( ); j++ ) {
	const basevector& e = local_HBVs[i].EdgeObject(j);
	if ( pass == 1 ) nbases += e.size( );
	else {
	  HKP_bases.push_back( e );
	  base_to_HKP_ID.push_back( i );
	}
      }
    }
    
    if ( pass == 1 ) {
      HKP_bases.Reserve( nbases/16 + n_HKP_paths, n_HKP_paths );
      base_to_HKP_ID .reserve( n_HKP_paths );
    }
  }
  
  // Add to this basevector all of the global unibases.  Now it contains the
  // entire dataset of basevectors known to be a part of the assembly.
  vecbasevector bases = HKP_bases;
  bases.Append( unibases );
  
  // Re-path.  This is a runtime bottleneck.
  cout << Date( ) << ": Pathing this vecbasevector" << endl;
  bool parallel = (NUM_THREADS > 1);
  vecKmerPath spaths, spaths_rc;
  vec<tagged_rpint> spathsdb;
  ReadsToPathsCoreY(bases, K, spaths, spaths_rc, spathsdb, 
                    run_dir + "/MergeNeighborhoods", NUM_THREADS);
  
  // Write the bases and paths to file.  EvalHyper needs this data, and it's
  // smarter to output it here than to force EvalHyper to re-path everything.
  cout << Date( ) << ": Writing bases/paths to sub_dir" << endl;
  bases.WriteAll( sub_dir + "/reads.fastb" );
  spaths.WriteAll( sub_dir + "/reads.paths" + kK );
  spaths_rc.WriteAll( sub_dir + "/reads.paths_rc" + kK );
  BinaryWriter::writeFile( sub_dir + "/reads.pathsdb" + kK, spathsdb );
  
  
  // Convert the local HyperBasevectors into HyperKmerPaths that use the
  // newly generated kmer-numbering system.
  cout << Date( ) << ": Converting local HyperBasevectors to HyperKmerPaths, and merging" << endl;
  int count = 0;
  vec<HyperKmerPath> local_HKPs;
  for ( int i = 0; i < n_seeds; i++ ) {
    vec<KmerPath> these_paths( 0 );
    for ( int j = 0; j < local_HBVs[i].EdgeObjectCount( ); j++ )
      these_paths.push_back( spaths[count++] );
    HyperKmerPath h( K, local_HBVs[i], these_paths );
    local_HKPs.push_back(h);
  }
  BinaryWriter::writeFile( sub_dir + "/nhood.hypers", local_HKPs );

  // Merge the local HyperKmerPaths into a single HyperKmerPath.
  HyperKmerPath hkp( K, local_HKPs );
  
  // Find "unique" unipaths and make a pathsdb out of them.
  vecKmerPath uniq_unipathsx;
  for ( size_t i = 0; i < n_unipaths; i++ ) {
    // Unipaths are marked as unique if their copy number is in the range
    // [1, ploidy] and if they are longer than 20 kmers.
    if ( predicted_CNs[i] <= ploidy && predicted_CNs[i] > 0 &&
	 unibases[i].size( ) - K + 1 >= 20 )
      uniq_unipathsx.push_back_reserve( spaths[ n_HKP_paths + i ] );
  }
  
  vec<tagged_rpint> uniqdb;
  CreateDatabase( uniq_unipathsx, uniqdb );
  
  // Merge the HyperKmerPath via GroupedInternalMerge.
  // This is the central algorithm of MergeNeighborhoods, and it is a serious
  // runtime bottleneck.
  KmerBaseBroker kbb( K, spaths, spaths_rc, spathsdb, bases );
  NegativeGapValidator ngv(&kbb);
  cout << Date( ) << ": Calling GroupedInternalMergeLG to merge " << n_seeds
       << " local neighborhoods" << endl;
  GroupedInternalMergeLG( local_HKPs, hkp, HKP_bases, base_to_HKP_ID, sub_dir,
			  NUM_THREADS,
			  MAX_KMER_FREQ,
			  MIN_ALIGN_LENGTH, MAX_GROUP_SIZE,
			  MIN_OVERLAP, MIN_PROPER_OVERLAP,
			  MIN_PROPER_OVERLAP_FINAL,
			  MAX_EXTENSION_QUEUE_SIZE,
			  &ngv, uniqdb, CHECKPOINT_INTERVAL,
			  sub_dir + "/hyper.ckpt", SKIP_FINAL_MERGE );
  
  cout << Date( ) << ": Cleaning up assembly HyperKmerPath" << endl;
  HKPCleanup(hkp);

  // Write assembly to file.
  cout << Date( ) << ": Writing assembly to file <sub_dir>/hyper" << endl;
  BinaryWriter::writeFile( sub_dir + "/hyper", hkp );
  Ofstream( dot, sub_dir + "/hyper.dot" );
  hkp.PrintSummaryDOT0w( dot, True, False, True, 0, True );
  
  cout << Date( ) << ": Done!" << endl;
  return 0;
}
