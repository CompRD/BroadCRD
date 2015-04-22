// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#include "assembly/CoverageAnnotator.h"


void CoverageAnnotator::Annotate( Contig& theContig, 
                                  vec<ContigAnnotation>& annotations,
                                  bool append ) const
{
  if ( ! append )
    annotations.clear();

  vec<ReadLocation> readLocs;
  theContig.GetReadLocations( readLocs );
  sort( readLocs.begin(), readLocs.end() );

  vec<int> begins, ends;
  for( unsigned int locIdx = 0; locIdx < readLocs.size(); ++locIdx )
  {
    begins.push_back( max( readLocs[locIdx].Begin(), 0 ) );
    ends.push_back( min( readLocs[locIdx].End(), theContig.GetLength() ) );
  }

  int coverage = 0;
  int regionStart = 0;

  unsigned int beginIdx = 0, endIdx = 0;

  while ( beginIdx < begins.size() || endIdx < ends.size() )
  {
    if ( beginIdx == begins.size() ||
         ends[endIdx] <= begins[beginIdx] )
    {
      if ( coverage == m_lowThres )
      {
        // Entering low coverage.
        regionStart = ends[ endIdx ];
      }

      --coverage;

      if ( coverage == m_highThres )
      {
        // Leaving high coverage.
        if ( ends[endIdx] - regionStart > 0 )
          annotations.push_back( ContigAnnotation( theContig, Annotation( "high coverage" ),
                                                   Interval( regionStart, ends[endIdx] ) ) );
        regionStart = ends[endIdx];
      }

      ++endIdx;
    }
    else
    {
      if ( coverage == m_highThres )
      {
        // Entering high coverage.
        regionStart = begins[ beginIdx ];
      }

      ++coverage;

      if ( coverage == m_lowThres )
      {
        // Leaving low coverage.
        if ( begins[beginIdx] - regionStart > 0 )
          annotations.push_back( ContigAnnotation( theContig, Annotation( "low coverage" ),
                                                   Interval( regionStart, begins[beginIdx] ) ) );
        regionStart = begins[beginIdx];
      }

      ++beginIdx;
    }
  }

  // We always need to check for a low coverage region at the end.
  if ( theContig.GetLength() - regionStart > 0 )
    annotations.push_back( ContigAnnotation( theContig, Annotation( "low coverage" ),
                                             Interval( regionStart, theContig.GetLength() ) ) );
}
