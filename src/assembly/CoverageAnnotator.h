// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#ifndef ASSEMBLY_COVERAGEANNOTATOR_H
#define ASSEMBLY_COVERAGEANNOTATOR_H

#include "assembly/AssemblyContig.h"
#include "assembly/ContigAnnotation.h"

class CoverageAnnotator
{
 public:

  CoverageAnnotator( const int lowThres = 0,
                     const int highThres = 10000000 )
    : m_lowThres( lowThres ),
      m_highThres( highThres )
  {}

  // Regions of contigs with coverage less than this will be annotated
  // as "low coverage".
  void SetLowThreshold( const int lowThres )
  {
    m_lowThres = lowThres;
    ForceAssertLt( m_lowThres, m_highThres );
  }

  // Regions of contigs with coverage greater than this will be annotated
  // as "high coverage".
  void SetHighThreshold( const int highThres )
  {
    m_highThres = highThres;
    ForceAssertLt( m_lowThres, m_highThres );
  }

  // Check the coverage of the given contig by its reads, creating
  // annotations for contiguous regions of "low coverage" and "high
  // coverage".  These will be appended in sorted order to the
  // annotations vector.
  void Annotate( Contig& theContig, 
                 vec<ContigAnnotation>& annotations,
                 bool append = false ) const;

 private:
  int m_lowThres;
  int m_highThres;
};

#endif
