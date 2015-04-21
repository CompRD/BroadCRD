/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "LocsHandler.h"
#include "SeqInterval.h"
#include "SeqIntervalUtil.h"
#include "Vec.h"
#include "math/Functions.h"
#include "math/NStatsTools.h"

/**
 * IntersectIntervalExternal
 */
vec<int> IntersectIntervalExternal( const seq_interval &interval,
				    const lhandler &locs )
{
  vec<int> result;

  int contig_id = interval.SeqId( );
  int floc = locs.FirstLoc( contig_id );
  if ( floc < 0 )
    return result;
  int interval_a = interval.Begin( );
  int interval_b = interval.End( );
  for (int ii=floc; ii<locs.Size( ); ii++) {
    const read_location &loc = locs[ii];
    if ( loc.Contig( ) != contig_id )
      break;
    int read_a = loc.StartOnContig( );
    if ( read_a >= interval_b )
      break;
    int read_b = 1 + loc.StopOnContig( );
    if ( IntervalOverlap( read_a, read_b, interval_a, interval_b ) > 0 )
      result.push_back( ii );
  }

  return result;
}

/**
 * IntersectIntervalInternal
 */
vec<int> IntersectIntervalInternal( const seq_interval &interval,
				    const lhandler &locs )
{
  vec<int> result;

  int contig_id = interval.SeqId( );
  int floc = locs.FirstLoc( contig_id );
  if ( floc < 0 )
    return result;
  int interval_a = interval.Begin( );
  int interval_b = interval.End( );
  for (int ii=floc; ii<locs.Size( ); ii++) {
    const read_location &loc = locs[ii];
    if ( loc.Contig( ) != contig_id )
      break;
    int read_a = loc.StartOnContig( );
    if ( read_a >= interval_b )
      break;
    int read_b = 1 + loc.StopOnContig( );
    int read_len = loc.LengthOfRead( );
    int overlap = IntervalOverlap( read_a, read_b, interval_a, interval_b );
    if ( overlap == read_len )
      result.push_back( ii );
  }

  return result;
}

/**
 * FlattenSeqIntervals
 */
void FlattenSeqIntervals( const vec<seq_interval>& to_flatten,
			  vec<seq_interval>& ans )
{
  ans.clear(); ans.reserve( to_flatten.size() );
  for( vec<seq_interval>::const_iterator i = to_flatten.begin();
       i != to_flatten.end(); i++ ) {
    if( ans.empty()
	|| i->SeqId() != ans.back().SeqId()
	|| i->Begin() > ans.back().End() )
      ans.push_back( *i );
    else if( i->End() > ans.back().End() )
      ans.back().SetEnd( i->End() );
  }
}

/**
 * SelectIntersecting
 */
void SelectIntersecting( unsigned long n_contigs,
			 const vec<seq_interval>& in_wins,
			 const vec<seq_interval>& annotations,
			 vec<seq_interval> &out_wins )
{
  ForceAssert( is_sorted( annotations.begin( ), annotations.end( ) ) );

  out_wins.clear( );
  out_wins.reserve( in_wins.size( ) );

  vec<int> fannots( n_contigs, -1 );
  for (int ii=annotations.size( )-1; ii>=0; ii--)
    fannots[ annotations[ii].SeqId( ) ] = ii;

  for (int ii=0; ii<(int)in_wins.size( ); ii++) {
    const seq_interval &win = in_wins[ii];
    int contig_id = win.SeqId( );
    int fannot = fannots[contig_id];
    if ( fannot < 0 )
      continue;
    for (int jj=fannot; jj<(int)annotations.size( ); jj++) {
      if ( annotations[jj].SeqId( ) != contig_id )
	break;
      if ( annotations[jj].HasOverlapWith( win ) ) {
	out_wins.push_back( win );
	break;
      }
    }
  }

}

/**
 * PrintSeqIntBasicNStats
 */
void PrintSeqIntBasicNStats( const vec<seq_interval> &intervals,
			     const String &name,
			     ostream &out )
{
  vec<int> lengths;
  lengths.reserve( intervals.size( ) );
  for (uint ii=0; ii<intervals.size( ); ii++) {
    int len = intervals[ii].Length( );
    if ( len > 0 ) lengths.push_back( len );
  }

  PrintBasicNStats( name, lengths, out );
}

