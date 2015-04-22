///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "paths/KmerPath.h"
#include "paths/KmerPathInterval.h"
#include "paths/assisted/KmerPathPerfectOverlaps.h"

/**
 * GlueKmerPaths
 *
 * It does not test that the two reads overlap, it just glues them at
 * the specified segments.
 */
KmerPath GlueKmerPaths( const KmerPath &left,
			const KmerPath &right,
			const int lpos,
			const int rpos )
{
  ForceAssert( lpos < left.NSegments( ) );
  ForceAssert( rpos < right.NSegments( ) );

  KmerPath answer;
  answer.reserve( lpos + right.NSegments( ) - rpos );
  
  for (int ii=0; ii<lpos; ii++)
    answer.AddSegment( left.Segment( ii ) );

  const KmerPathInterval &lseg = left.Segment( lpos );
  const KmerPathInterval &rseg = right.Segment( rpos );
  longlong jstart = Min( lseg.Start( ), rseg.Start( ) );
  longlong jstop = Max( lseg.Stop( ), rseg.Stop( ) );
  answer.AddSegment( jstart, jstop );

  for (int ii=rpos+1; ii<right.NSegments( ); ii++)
    answer.AddSegment( right.Segment( ii ) );

  return answer;
}

/**
 * PerfectOverlap
 */
bool PerfectOverlap( const KmerPath &rpath,
		     const KmerPath &upath,
		     const int rpos,
		     const int upos )
{
  int r_segs = rpath.NSegments( );
  int u_segs = upath.NSegments( );

  // Look for the beginning of the align.
  int rseg = rpos;
  int useg = upos;
  while ( rseg > 0 && useg > 0 ) {
    rseg--;
    useg--;
  }
  
  // Loop over all segments.
  while ( useg < u_segs && rseg < r_segs ) {
    const KmerPathInterval &r_kpi = rpath.Segment( rseg );
    const KmerPathInterval &u_kpi = upath.Segment( useg );
    if ( ! r_kpi.Overlaps( u_kpi ) ) return false;

    // Left.
    if ( rseg == 0 && useg > 0 )
      if ( r_kpi.Start( ) < u_kpi.Start( ) ) return false;
    if ( useg == 0 && rseg > 0 )
      if ( u_kpi.Start( ) < r_kpi.Start( ) ) return false;
    if ( rseg > 0 && useg > 0 )
      if ( r_kpi.Start( ) != u_kpi.Start( ) ) return false;
    
    // Right.
    if ( rseg == r_segs - 1 && useg < u_segs - 1 )
      if ( r_kpi.Stop( ) > u_kpi.Stop( ) ) return false;
    if ( useg == u_segs - 1 && rseg < r_segs - 1 )
      if ( u_kpi.Stop( ) > r_kpi.Stop( ) ) return false;
    if ( rseg < r_segs - 1 && useg < u_segs - 1 )
      if ( r_kpi.Stop( ) != u_kpi.Stop( ) ) return false;

    useg++;
    rseg++;
  }
  
  // Is perfect overlap.
  return true;
  
}

/**
 * IsReadEmbedded
 */
bool IsReadEmbedded( const KmerPath &rpath,
		     const KmerPath &upath,
		     const int rpos,
		     const int upos )
{
  // Must be a perfect overlap.
  if ( ! PerfectOverlap( rpath, upath, rpos, upos ) ) return false;

  // Look for the beginning of the align.
  int r_segs = rpath.NSegments( );
  int u_segs = upath.NSegments( );
  int rseg = rpos;
  int useg = upos;
  while ( rseg > 0 && useg > 0 ) {
    rseg--;
    useg--;
  }

  // First segment mismatch.
  if ( rseg != 0 ) return false;

  const KmerPathInterval *r_kpi = &( rpath.Segment( 0 ) );
  const KmerPathInterval *u_kpi = &( upath.Segment( useg ) );
  if ( r_kpi->Start( ) < u_kpi->Start( ) ) return false;

  // Go to end of align.
  rseg = rpos;
  useg = upos;
  while ( rseg < r_segs-1 && useg < u_segs-1 ) {
    rseg++;
    useg++;
  }

  // Last segment mismatch.
  if ( rseg != r_segs - 1 ) return false;

  r_kpi = &( rpath.Segment( rseg ) );
  u_kpi = &( upath.Segment( useg ) );
  if ( r_kpi->Stop( ) > u_kpi->Stop( ) ) return false;
  
  // Is embedded.
  return true;
  
}

/**
 * IsProperExtension
 */
bool IsProperExtension( const KmerPath &rpath,
			const KmerPath &upath,
			const int rpos,
			const int upos )
{
  // Must be a perfect overlap.
  if ( ! PerfectOverlap( rpath, upath, rpos, upos ) ) return false;

  // Look for the end of the align.
  int r_segs = rpath.NSegments( );
  int u_segs = upath.NSegments( );
  int rseg = rpos;
  int useg = upos;
  while ( rseg < r_segs-1 && useg < u_segs-1 ) {
    rseg++;
    useg++;
  }
  
  // Last segment.
  const KmerPathInterval & r_lastseg = rpath.Segment( rseg );
  const KmerPathInterval & u_lastseg = upath.Segment( useg );
  if ( useg != u_segs-1 ) return false;
  if ( rseg == r_segs-1 ) return ( u_lastseg.Stop( ) < r_lastseg.Stop( ) );
  
  // Is proper.
  return true;

}

/**
 * SimplePrintAlign
 */
void SimplePrintAlign( const KmerPath &path1,
		       const KmerPath &path2,
		       const int pos1,
		       const int pos2,
		       ostream &out )
{
  // Look for the beginning of the align.
  int nsegs1 = path1.NSegments( );
  int nsegs2 = path2.NSegments( );
  int seg1 = pos1;
  int seg2 = pos2;
  while ( seg1 > 0 && seg2 > 0 ) {
    seg1--;
    seg2--;
  }
  
  // Build printable table.
  vec< vec<String> > table( 2 );
  while ( seg1 < nsegs1 && seg2 < nsegs2 ) {
    String openb =  ( seg1 == pos1 ) ? "{" : "[";
    String closeb =  ( seg1 == pos1 ) ? "}" : "]";

    // Segment on first kmer path.
    if ( 0 <= seg1 && seg1 < nsegs1 ) {
      const KmerPathInterval &kpi = path1.Segment( seg1 );
      table[0].push_back( openb
			  + ToString( kpi.Start( ) )
			  + "-"
			  + ToString( kpi.Stop( ) )
			  + closeb );
    }
    else
      table[0].push_back( "" );

    // Segment on second kmer path.
    if ( 0 <= seg2 && seg2 < nsegs2 ) {
      const KmerPathInterval &kpi = path2.Segment( seg2 );
      table[1].push_back( openb
			  + ToString( kpi.Start( ) )
			  + "-"
			  + ToString( kpi.Stop( ) )
			  + closeb );
    }
    else table[1].push_back( "" );
    
    // Next segment.
    seg1++;
    seg2++;
  }
  
  // Print table.
  PrintTabular( cout, table, 3 );
  
}
