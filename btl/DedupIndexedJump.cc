///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"

#include "Fastavector.h"
#include "Intvector.h"
#include "TokenizeString.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"
#include "pairwise_aligners/KmerAligner.h"
#include "pairwise_aligners/SmithWatBandedA.h"

#include <omp.h>
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

/**
 * DedupCore
 *
 * It returns the number of duplicated (and removed) reads.
 */
int DedupCore( const vecbvec &reads,
	       const IntVec& ids,
	       KmerAligner<24> &aligner24,
	       IntVec &select )
{
  // HEURISTICS.
  const int band = 2;
  const int max_mis = 5;
  const int max_ind = 2;
  
  // Clean up.
  select.clear( );
  int dups_found = 0;

  // This should not happen.
  if ( ids.size( ) < 1 ) return 0;
  
  // Build select.
  select.push_back( ids[0] );
  for (size_t ii=1; ii<ids.size( ); ii++) {
    vec<int> kapos;
    vec<int> kaids;
    aligner24.FindPossibleAlignments( reads[ ids[ii] ], kapos, kaids );

    bool is_dup = false;
    for (size_t jj=0; jj<kaids.size( ); jj++) {
      IntVec::iterator it = find( select.begin( ), select. end( ), kaids[jj] );
      if ( it == select.end( ) ) continue;
      
      // reads[ ids[ii] ] may be a duplicate of reads[*it].
      const int off = - kapos[jj];
      
      align al;
      int err = 0;
      SmithWatBandedA( reads[ ids[ii] ], reads[*it], off, band, al, err );
      if ( al.Pos1( ) - al.pos1( ) < 24 ) continue;
      
      unsigned int iilen = reads[ ids[ii] ].size( );
      unsigned int itlen = reads[*it].size( );
      int n_err = al.Errors( reads[ ids[ii] ], reads[*it] );

      look_align theLA( ids[ii], *it, iilen, itlen, false, al, 0, n_err, 0);
      if ( ! theLA.IsProper( ) ) continue;
      
      int n_muts = theLA.Mutations( );
      int n_ind = theLA.Indels( );
      bool good = ( n_muts <= max_mis && n_ind <= max_ind );
      if ( n_muts > max_mis || n_ind > max_ind ) continue;
      
      // A duplicate!
      is_dup = true;
      dups_found++;
      break;
    }
    if ( is_dup ) continue;
    
    // Add ids[ii] to select.
    select.push_back( ids[ii] );
    
  }
  
  // Done.
  return dups_found;

}

/**
 * DedupIndexedJump
 *
 * It will find and remove duplicate reads with the same index.
 *
 * FASTA: full path name to reads' fasta file.
 * MAP: full path name to map file (on each line: "read_name\tbarcode_id").
 * OUT: full path name of output file (same format as MAP).
 * MAX_BAG_SIZE: skip bags with more than these many reads.
 * NUM_THREADS: use all if 0.
 *
 * REMARK: this module is being replaced by Dedup( ) in CBarcodes.h,
 *         it is only kept for backward compatibility with some of
 *         Aaron's scripts (for testing purposes).
 */
int main( int argc, char *argv[] )
{
  RunTime( );
 
  BeginCommandArguments;
  CommandArgument_String( FASTA );
  CommandArgument_String( MAP );
  CommandArgument_String( OUT );
  CommandArgument_Int_OrDefault( MAX_BAG_SIZE, 1000 );
  CommandArgument_UnsignedInt_OrDefault( NUM_THREADS, 0 );
  EndCommandArguments;
  
  // Thread control.
  NUM_THREADS = configNumThreads( NUM_THREADS );
  omp_set_num_threads( NUM_THREADS );
  
  // Load.
  cout << Date( ) << ": loading reads... " << flush;
  vec<fastavector> readsfa;
  vec<String> tmpnames;
  LoadFromFastaFile( FASTA, readsfa, tmpnames );
  const size_t n_reads = readsfa.size( );

  vecString names;
  names.reserve( n_reads );
  for (size_t ii=0; ii<n_reads; ii++)
    names.push_back( tmpnames[ii] );

  vecbvec reads;
  reads.reserve( n_reads );
  for (size_t ii=0; ii<n_reads; ii++)
    reads.push_back( readsfa[ii].ToBasevector( ) );

  IntVec ids;
  ids.reserve( n_reads );
  for (size_t ii=0; ii<n_reads; ii++)
    ids.push_back( (longlong)ii );

  vecString names_orig = names;
  SortSync( names, ids );

  vec< pair<String,longlong> > bar2ids;
  bar2ids.reserve( reads.size( ) );
  ifstream map_in( MAP.c_str( ) );
  while ( map_in ) {
    String line;
    vec<String> tokens;
    getline( map_in, line );
    if ( ! map_in ) break;
    Tokenize( line, '\t', tokens );
    ForceAssertEq( (int)tokens.size( ), 2 );

    vecString::iterator it = lower_bound(names.begin(), names.end(), tokens[0]);
    ForceAssert( it != names.end( ) );
    ForceAssertEq( *it, tokens[0] );
    longlong id = distance( names.begin( ), it );

    bar2ids.push_back( make_pair( tokens[1], ids[id] ) );
  }
  map_in.close( );
  sort( bar2ids.begin( ), bar2ids.end( ) );
  size_t n_bars = bar2ids.size( );
  cout << ToStringAddCommas( n_reads ) << " found" << endl;

  // No mapping info?
  if ( n_bars < 1 ) {
    cout << "\nFatal errors: no mapping info found.\n" << endl;
    return 1;
  }

  // Create bags of reads.
  cout << Date( ) << ": building bags... " << flush;
  size_t n_bags = 1;
  for (size_t ii=1; ii<n_bars; ii++)
    if ( bar2ids[ii].first != bar2ids[ii-1].first )
      n_bags++;

  vecString bag_names;
  VecIntVec bag_rids;
  bag_names.reserve( n_bags );
  bag_rids.reserve( n_bags );

  for (size_t ii=0; ii<n_bars; ii++) {
    const String &bag_name = bar2ids[ii].first;
    const longlong rid = bar2ids[ii].second;
    if ( ii == 0 || bag_name != bag_names.back( ) ) {
      bag_names.push_back( bag_name );
      bag_rids.push_back( IntVec( ) );
    }
    bag_rids[ bag_rids.size( )-1 ].push_back( rid );
  }
  cout << ToStringAddCommas( n_bags ) << " found" << endl;

  // Kmerizing reads.
  cout << Date( ) << ": kmerizing reads (using k=24)" << endl;
  KmerAligner<24> aligner24;
  aligner24.SetBases( reads );

  // Output stream.
  ofstream out( OUT.c_str( ) );

  // Log progress.
  const size_t n_dotter = 10000;
  cout << Date( ) << ": digesting "
       << ToStringAddCommas( n_bags ) << " bags"
       << endl;

  // Loop over all bags.
  int totd = 0;
  int nhuge = 0;
  int ndone = 0;
  #pragma omp parallel for
  for( size_t bag_id=0; bag_id<n_bags; bag_id++) {
    
    if ( (int)bag_rids[bag_id].size( ) > MAX_BAG_SIZE ) {
      #pragma omp critical
      {
	nhuge++;
	ndone++;
	for (size_t ii=0; ii<bag_rids[bag_id].size( ); ii++)
	  out << names_orig[ bag_rids[bag_id][ii] ] << "\t"
	      << bag_names[bag_id] << "\n";
      }
      continue;
    }

    #pragma omp critical
    if ( ndone > 0 && ndone % n_dotter == 0 ) {
      double ratio = 100. * SafeQuotient( ndone, n_bags );
      cout << "  "
	   << Date( ) << ": "
	   << ToStringAddCommas( ndone ) << " bags done ("
	   << ToString( ratio, 1 ) << "%) "
	   << totd << " dups found"
	   << endl;
    }

    IntVec select;
    int ndups = DedupCore( reads, bag_rids[bag_id], aligner24, select );
    
    #pragma omp critical
    {
      totd += ndups;
      ndone++;
      for (size_t ii=0; ii<select.size( ); ii++)
	out << names_orig[ select[ii] ] << "\t" << bag_names[bag_id] << "\n";
    }
  }

  // Print warning, if huge bag were found.
  if ( nhuge > 0 )
    cout << "\nWARNING: " << ToStringAddCommas( nhuge )
	 << " bags had more than " << ToStringAddCommas( MAX_BAG_SIZE )
	 << " reads, and were not tested\n" << endl;

  // Close stream (and report).
  cout << Date( ) << ": " << ToStringAddCommas( totd ) << " dups found" << endl;
  out.close( );
  
  // Done.
  cout << Date( ) << ": DedupIndexedJump done" << endl;

}
