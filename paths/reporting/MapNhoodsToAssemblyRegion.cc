///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "lookup/LookAlign.h"
#include "paths/ReadHBVs.h"
#include "paths/reporting/MapNhoodsUtils.h"
#include "util/RunCommand.h"

#include <omp.h>
// MakeDepend: library OMP
// MakeDepend: dependency MakeLookupTable

/**
 * MapNhoodsToAssemblyRegion
 *
 * Interactive tool to map nhoods onto a region of a given assembly.
 * To do so, it first generates a rough estimate of where nhoods land
 * (by aligning the seeds to the whole assembly), and then it
 * realignes the candidate nhoods to the local region.
 *
 * Interactive mode: select ranges of the type SUPER_ID:BEGIN-END
 *
 * Non interactive modes:
 *   a. specify a range (in the same format as in the interactive mode)
 *   b. specify one ore more nhood ids
 *
 * Warning: the first time the code is run, it will generate and cache
 * several files - this is a time consuming part.
 *
 * READS: it loads ../<READS>.unibases.k<K>
 * ASSEMBLY: assembly head
 * OUT_DIR: relative to <SUBDIR>
 * RANGE: if given, just do this range and exit (not interactive mode)
 * NHOOD_IDS: if given, just align these nhoods and exit (not interactive mode)
 * SKIP_UNMAPPED: skip nhoods with unmapped seeds
 * FORCE: do not used cached data
 * NUM_THREADS: use all available if 0
 * ORIGIN: shift starts on target by this amount
 */
int main( int argc, char *argv[] )
{
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_Int( K );
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String_OrDefault( READS, "all_reads" );
  CommandArgument_String_OrDefault( ASSEMBLY, "final" );
  CommandArgument_String_OrDefault( OUT_DIR, ASSEMBLY + ".nhoods.map" );
  CommandArgument_String_OrDefault( RANGE, "" );
  CommandArgument_String_OrDefault( NHOOD_IDS, "" );
  CommandArgument_Bool_OrDefault( SKIP_UNMAPPED, True );
  CommandArgument_Bool_OrDefault( FORCE, False );
  CommandArgument_Int_OrDefault( NUM_THREADS, 0 );
  CommandArgument_Int_OrDefault( ORIGIN, 0 );
  EndCommandArguments;

  // Dir and file names.
  String data_dir = PRE + "/" + DATA;
  String run_dir = data_dir + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  String out_dir = sub_dir + "/" + OUT_DIR;
  
  String strK = ToString( K );
  String unibases_file = run_dir + "/" + READS + ".unibases.k" + strK;
  String seeds_ids_file = sub_dir + "/seeds.ids";
  String assembly_head = sub_dir + "/" + ASSEMBLY;
  String assembly_lookup_head = assembly_head + ".assembly";
  String assembly_lookup_file = assembly_head + ".assembly.lookup";
  String assembly_fasta_file = assembly_head + ".assembly.fasta";
  String assembly_lengths_file = assembly_head + ".assembly.lengths";
  
  String seeds_aligns_file = out_dir + "/seeds.qlt";
  String seeds_bases_file = out_dir + "/seeds_bases.fastb";
  
  Mkpath( out_dir );
  
  // Thread control
  NUM_THREADS = configNumThreads( NUM_THREADS );
  omp_set_num_threads( NUM_THREADS );

  // Generate lookup table.
  if ( FORCE || ! IsRegularFile( assembly_lookup_file ) ) {
    cout << Date( ) << ": generating lookup table" << endl;
    String theCommand
      = "MakeLookupTable SOURCE=" + assembly_fasta_file
      + " OUT_HEAD=" + assembly_lookup_head;
    RunCommandWithLog( theCommand, "/dev/null" );
  }
  
  // Seed ids must be sorted.
  READ( seeds_ids_file, vec<int>, seeds_ids );
  ForceAssert( is_sorted( seeds_ids.begin( ), seeds_ids.end( ) ) );
  int n_seeds = seeds_ids.size( );
  
  // Align seeds.
  vecbvec seeds_bases;
  if ( FORCE || ! IsRegularFile( seeds_aligns_file ) ) {
    cout << Date( ) << ": loading " << n_seeds << " seeds" << endl;
    seeds_bases.SparseRead( unibases_file, seeds_ids, 0 );
    seeds_bases.WriteAll( seeds_bases_file );
  }
  
  cout << Date( ) << ": aligning seeds... " << flush;
  vec<look_align> seeds_aligns;
  GetAlignsFast( K, seeds_bases_file, assembly_lookup_file, seeds_aligns_file,
		 seeds_aligns, ! FORCE, out_dir );
  cout << seeds_aligns.size( ) << " aligns found" << endl;

  // Build an index for the aligned seeds (-1: unmapped, -2: multiplet).
  vec<int> aligns_index( n_seeds, -1 );
  for (int ii=0; ii<seeds_aligns.isize( ); ii++) {
    int seed_id = seeds_aligns[ii].query_id;
    vec<int>::iterator it
      = lower_bound( seeds_ids.begin( ), seeds_ids.end( ), seed_id );
    int nhood_id = distance( seeds_ids.begin( ), it );
    ForceAssertLt( nhood_id, seeds_ids.isize( ) );
    ForceAssertEq( seeds_ids[nhood_id], seed_id );
    if ( aligns_index[nhood_id] == - 1 ) aligns_index[nhood_id] = ii;
    else aligns_index[nhood_id] = -2;
  }
  
  // Load assembly.
  cout << Date( ) << ": loading global assembly" << endl;
  vec<int> assembly_lengths;
  vec<fastavector> assembly_fasta;
  LoadFromFastaFile( assembly_fasta_file, assembly_fasta );
  assembly_lengths.reserve( assembly_fasta.size( ) );
  for (int ii=0; ii<assembly_fasta.isize( ); ii++)
    assembly_lengths.push_back( assembly_fasta[ii].size( ) );
  
  // Load local assemblies.
  cout << Date( ) << ": loading local assemblies" << endl;
  vec<bool> all_seeds( seeds_ids.size( ), true );
  vec<HyperBasevector> hypers;
  readHBVs( sub_dir, all_seeds, &hypers, false );
  
  cout << endl;

  // Non interactive mode - by RANGE.
  if ( RANGE != "" ) {
    MapRange( cout, RANGE, SKIP_UNMAPPED, FORCE, K, n_seeds, out_dir,
	      assembly_head, assembly_fasta, assembly_lengths, seeds_ids,
	      aligns_index, seeds_aligns, hypers );
    
    cout << Date( ) << ": done" << endl;
    return 1;
  }
  
  // Non interactive mode - by NHOOD_IDS.
  if ( NHOOD_IDS != "" ) {
    vec<int> select;
    ParseIntSet( NHOOD_IDS, select );
    for (int ii=0; ii<select.isize( ); ii++) {
      cout << ii + 1 << "." << select.size( ) << flush;

      int nhood_id = select[ii];
      triple<int,int,int> placement( -1, -1, -1 );
      const HyperBasevector &hyper = hypers[nhood_id];
      const look_align &al = seeds_aligns[ aligns_index[nhood_id] ];
      placement = EstimatePlacement( hyper, al, assembly_lengths );

      vec<look_align> aligns;
      RefinePlacement( placement, aligns, K, nhood_id, assembly_head,
		       out_dir, assembly_fasta, hypers, FORCE );

      cout << "\tnhood " << nhood_id
	   << "\tu_" << seeds_ids[nhood_id]
	   << "\ton t_" << placement.first
	   << " [" << placement.second
	   << ", " << placement.third
	   << ")_" << assembly_lengths[placement.first]
	   << endl;
    }

    cout << "\n" << Date( ) << ": done" << endl;
    return 1;      
  }
  
  // Call MapRange (interactive mode).
  while ( 1 ) {
    cout << "> " << flush;
    
    String cmnd;
    getline( cin, cmnd );
    if ( ! cin || cmnd == "q" ) {
      cout << endl;
      break;
    }
    
    MapRange( cout, cmnd, SKIP_UNMAPPED, FORCE, K, n_seeds, out_dir,
	      assembly_head, assembly_fasta, assembly_lengths, seeds_ids,
	      aligns_index, seeds_aligns, hypers );
  }
  
  // Done.
  cout << Date( ) << ": done" << endl;

}
