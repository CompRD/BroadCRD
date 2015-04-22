///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Vec.h"
#include "btl/ClusterOffsets.h"
#include "math/Functions.h"

/**
 * ClusterOffsets
 */
void ClusterOffsets( const int MIN_SEEDS,
		     const vec<int> &offsets,
		     vec<int> &clusters )
{
  ForceAssert( MIN_SEEDS > 0 );
  ForceAssert( is_sorted( offsets.begin( ), offsets.end( ) ) );
  clusters.clear( );

  // HEURISTICS.
  const float MIN_MEDIAN = 10.0;
  const float MEDIAN_MULTIPLIER = 5.0;

  // Nothing to do.
  if ( offsets.isize( ) < MIN_SEEDS ) return;
  
  // Special case: offsets has one entry.
  if ( offsets.size( ) == 1 ) {
    clusters.resize( 1, offsets[0] );
    return;
  }

  // Gather distances between offsets.
  vec<float> deltas;
  deltas.reserve( offsets.size( ) - 1 );
  for (int ii=1; ii<offsets.isize( ); ii++)
    deltas.push_back( offsets[ii] - offsets[ii-1] );

  // No deltas found, leave.
  if ( deltas.size( ) < 1 ) return;

  // Median of deltas, define big gap size.
  sort( deltas.begin( ), deltas.end( ) );
  float med = Max( MIN_MEDIAN, Median( deltas ) );
  int max_gap = med * MEDIAN_MULTIPLIER;
  
  // Build clusters.
  vec<int> newc( 1, offsets[0] );
  for (int ii=1; ii<offsets.isize( ); ii++) {
    if ( offsets[ii] - offsets[ii-1] > max_gap ) {
      if ( newc.isize( ) >= MIN_SEEDS ) {
	sort( newc.begin( ), newc.end( ) );
	clusters.push_back( Median( newc ) );
      }
      newc.clear( );
    }
    newc.push_back( offsets[ii] );
  }
  if ( newc.isize( ) >= MIN_SEEDS ) {
    sort( newc.begin( ), newc.end( ) );
    clusters.push_back( Median( newc ) );    
  }
  
}
