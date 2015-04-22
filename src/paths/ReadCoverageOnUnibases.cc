///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "math/Functions.h"
#include "math/NStatsTools.h"
#include "paths/KmerPath.h"
#include "paths/Unipath.h"
#include "util/RunCommand.h"

/**
 * ReadCoverageOnUnibases
 *
 * Compute mean and stdev of read coverage of unipaths. Output is
 * saved in <HEAD>.unipaths.readcov.k<K> (a vec<normal_distribution>).
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_Int( K );
  CommandArgument_String( HEAD );
  EndCommandArguments;

  // File names.
  String strK = ToString( K );
  String pathsdb_file = HEAD + ".pathsdb.k" + strK;
  String unipaths_file = HEAD + ".unipaths.k" + strK;
  String readcov_file = HEAD + ".unipaths.readcov.k" + strK;
  
  // Check files exist.
  vec<String> needed;
  needed.push_back( pathsdb_file );
  needed.push_back( unipaths_file );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;

  // Load.
  cout << Date( ) << ": loading pathsdb" << endl;
  BREAD2( pathsdb_file, vec<tagged_rpint>, pathsdb );

  cout << Date( ) << ": loading unipaths" << endl;
  vecKmerPath unipaths( unipaths_file );

  cout << Date( ) << ": done loading" << endl;

  // Nstats of edge sizes.
  vec<int> lens;
  lens.reserve( unipaths.size( ) );
  for (size_t edge_id=0; edge_id<unipaths.size( ); edge_id++)
    lens.push_back( unipaths[edge_id].TotalLength( ) );
  cout << "\nEDGES N-STATS (KMER LENGTHS):\n\n";
  PrintBasicNStats( "", lens, cout );
  cout << endl;
  
  // Container for read coverages of edges (the answer).
  vec<NormalDistribution> read_covs;
  read_covs.reserve( unipaths.size( ) );

  // Loop over all edges.
  for (size_t edge_id=0; edge_id<unipaths.size( ); edge_id++) {
    const KmerPath &edge = unipaths[edge_id];

    // Container for coverage.
    vec<float> cov;
    cov.reserve( edge.TotalLength( ) );

    // Loop over all segments in edge.
    for (longlong segment_id=0; segment_id<edge.NSegments( ); segment_id++) {
      const KmerPathInterval &segment = edge.Segment( segment_id );
      const size_t start = segment.Start( );
      const size_t stop = segment.Stop( );

      // Loop over kmers in segment.
      for (size_t kmer_id=start; kmer_id<=stop; kmer_id++) {
	vec<longlong> locs;
	Contains( pathsdb, kmer_id, locs );
	cov.push_back( float( locs.size( ) ) );
      }
    }    
    
    // Add mean/stdev to answer.
    read_covs.push_back( SafeMeanStdev( cov ) );
  }
  
  // Save.
  cout << Date( ) << ": saving" << endl;
  WRITE( readcov_file, read_covs );

  // Done.
  cout << Date( ) << ": done" << endl;

}
