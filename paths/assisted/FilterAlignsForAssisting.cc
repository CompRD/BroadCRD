///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Vec.h"
#include "lookup/LookAlign.h"
#include "paths/assisted/FilterAlignsForAssisting.h"

/**
 * AlignsGap
 */
int AlignsGap( const look_align &left, const look_align &right )
{
  if ( left.target_id != right.target_id )
    return numeric_limits<int>::quiet_NaN( );
  
  int right_beg = right.pos2( ) - right.pos1( );
  int left_end = left.Pos2( ) + ( left.query_length - left.Pos1( ) );
  return right_beg - left_end;
}

/**
 * ChainsGap
 */
int ChainsGap( const vec<int> &left_chain,
	       const vec<int> &right_chain,
	       const vec<look_align> &aligns )
{
  const look_align &left = aligns[left_chain.back( )];
  const look_align &right = aligns[right_chain.front( )];

  return AlignsGap( left, right );
}

/**
 * FilterBigGaps
 */
void FilterBigGaps( const int MIN_GAP,
		    const int MAX_GAP,
		    vec<look_align> &aligns,
		    ostream *log )
{
  if ( log ) *log << "FILTER BY GAPS\n\n";
  
  // Sort aligns.
  order_lookalign_TargetBeginEnd sorter;
  if ( ! is_sorted( aligns.begin( ), aligns.end( ), sorter ) )
    sort( aligns.begin( ), aligns.end( ), sorter );
  
  // Build chains.
  vec< vec<int> > chains;
  for (int ii=0; ii<aligns.isize( ); ii++) {
    int gap = 0;
    if ( ii > 0 ) gap = AlignsGap( aligns[ii-1], aligns[ii] );

    bool is_new = ( chains.size( ) < 1 || gap < MIN_GAP || gap > MAX_GAP );
      
    if ( is_new ) chains.push_back( vec<int>( 1, ii ) );
    else chains[chains.size( )-1].push_back( ii );
  }

  // Tag chains for deletion.
  int n_deleted = 0;
  vec<bool> deleted( chains.size( ), false );
  for (int ii=1; ii<chains.isize( )-1; ii++) {
    int gap = ChainsGap( chains[ii-1], chains[ii+1], aligns );
    if ( MIN_GAP <= gap && gap <= MAX_GAP ) {
      deleted[ii] = true;
      n_deleted++;
      
      // Report event.
      if ( ! log ) continue;

      longlong n_bases_deleted = 0;
      for (int jj=0; jj<chains[ii].isize( ); jj++)
	n_bases_deleted += aligns[ chains[ii][jj] ].query_length;

      *log << "removing chain " << ii 
	   << " (" << chains[ii].size( ) << " align(s), "
	   << ToString( n_bases_deleted ) << " bases)\n";
    }
  }

  if ( log ) {
    if ( n_deleted < 1 ) *log << "no chain removed\n";
    *log << endl;
  }

  // Remove deleted aligns.
  vec<look_align> result;
  result.reserve( aligns.size( ) );
  for (int ii=0; ii<chains.isize( ); ii++) {
    if ( deleted[ii] ) continue;
    for (int jj=0; jj<chains[ii].isize( ); jj++)
      result.push_back( aligns[ chains[ii][jj] ] );
  }

  // Done.
  swap( aligns, result );

}

