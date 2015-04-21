///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "PairsManager.h"
#include "Vec.h"
#include "VecUtilities.h"
#include "paths/TranslateAligns.h"
#include "paths/reporting/PhysicalCoverageFromAlignsCore.h"

/**
 * PhysicalCoverageFromAlignsCore
 */
double PhysicalCoverageFromAlignsCore( ostream &log,
				       vec<alignlet> &aligns,
				       const vec<int> &index,
				       const vec<superb> &supers,
				       const PairsManager &pairs,
				       const double MAX_STRETCH,
				       const int RADIUS,
				       const bool VERBOSE )
{  
  // Convert aligns in super coordinates.
  TranslateAlignsToSupers( supers, aligns );
  
  // Core data (key: super_id).
  vec<double> total_insertlen( supers.size( ), .0 );
  vec<double> adjusted_superlen( supers.size( ), .0 );

  // The denumerators.
  double four_rad = 4.0 * double( RADIUS );
  for (int super_id=0; super_id<supers.isize( ); super_id++) {
    double slen = supers[super_id].FullLength( ) - four_rad;
    adjusted_superlen[super_id] = Max( .0, slen );
  }
  
  // The numerators (loop over all pairs).
  size_t n_valid_pairs = 0;
  for (size_t pair_id=0; pair_id<pairs.nPairs( ); pair_id++) {
    int id1 = pairs.ID1( pair_id );
    int id2 = pairs.ID2( pair_id );

    if ( index[id1] < 0 || index[id2] < 0 ) continue;
    
    const alignlet &al1 = aligns[ index[id1] ];
    const alignlet &al2 = aligns[ index[id2] ];
    if ( al1.TargetId( ) != al2.TargetId( ) ) continue;

    if ( al1.Fw1( ) == al2.Fw1( ) ) continue;

    const alignlet &alF = al1.Fw1( ) ? al1 : al2;
    const alignlet &alR = al1.Fw1( ) ? al2 : al1;
    int super_id = al1.TargetId( );
    int super_len = supers[super_id].FullLength( );
    if ( alF.pos2( ) < RADIUS ) continue;
    if ( alR.Pos2( ) > super_len - RADIUS ) continue;
    int begin = Max( alF.pos2( ), 2 * RADIUS );
    int end = Min( alR.Pos2( ), super_len - ( 2 * RADIUS ) );
    if ( end <= begin ) continue;
    
    total_insertlen[super_id] += double( end - begin );
    n_valid_pairs++;
  }

  // Find the coverages.
  vec<double> cov( supers.size( ), .0 );
  for (int super_id=0; super_id<supers.isize( ); super_id++) {
    if ( adjusted_superlen[super_id] < 1.0 ) continue;
    cov[super_id] = total_insertlen[super_id] / adjusted_superlen[super_id];
  }
  
  // Weighted mean.
  double total_slen = .0;
  for (int super_id=0; super_id<supers.isize( ); super_id++) {
    double slen = adjusted_superlen[super_id];
    if ( slen > .0 ) total_slen += slen;
  }
  
  if ( total_slen < 1.0 ) {
    cout << "FATAL ERROR: no supers left\n" << endl;
    return 0;
  }
      
  double mean_cov = .0;
  for (int super_id=0; super_id<supers.isize( ); super_id++) {
    if ( cov[super_id] < 1.0 ) continue;
    mean_cov += cov[super_id] * adjusted_superlen[super_id];
  }
  mean_cov = mean_cov / total_slen;
  
  // Print result.
  if ( VERBOSE ) {
    vec< vec<String> > table;
    vec<String> line;
    
    {
      line.push_back( "super_id" );
      line.push_back( "super_length" );
      line.push_back( "coverage" );
      table.push_back( line );
    }
    
    for (int super_id=0; super_id<supers.isize( ); super_id++) {
      line.clear( );
      line.push_back( ToString( super_id ) );
      line.push_back( ToString( supers[super_id].FullLength( ) ) );
      line.push_back( ToString( cov[super_id], 1 ) );
      table.push_back( line );
    }
    
    PrintTabular( log, table, 3, "rrr" );
  }
  
  log << "\n"
      << "SUMMARY\n"
      << "\n"
      << "pairs in input:   " << pairs.nPairs( ) << "\n"
      << "pairs aligned:    " << n_valid_pairs << "\n"
      << "overall coverage: " << ToString( mean_cov, 3 ) << "\n"
      << endl;

  // Return.
  return mean_cov;

}

