/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef C_SEEDS_HANDLER_H
#define C_SEEDS_HANDLER_H

#include "ReadLocation.h"
#include "String.h"
#include "STLExtensions.h"
#include "pairwise_aligners/CSeed.h"

/**
 * Class CSeedsHandler
 *
 * Organizer for seeds. As soon as all seeds are added, you must run
 * ClusterSeeds( ) to sort and generate the necessary maps.
 */
class CSeedsHandler {
  
public:
  
  CSeedsHandler( int read_id = -1 );

  void Add( int target_id, int start, Bool rc );

  // Organize seeds (must run this before any other querying method).
  void ClusterSeeds( );
  
  // Select cluster with more seeds.
  bool HeavySeed( CSeed &heavy_seed ) const;
  
  // Return the average of the seeds in the given cluster.
  CSeed AverageSeed( int fseedpos ) const;
  
  // One liner info.
  void PrintInfo( ostream &out ) const;
  
  
private:
  
  void SetDefaults( );
  
  
private:
  
  int read_id_;        // id of read
  int max_diff_;       // used to partition seeds into consistent bags
  vec<CSeed> seeds_;   // set of seeds
  vec<int> fseeds_;    // id in seeds_ of first seed in bag

};

#endif
