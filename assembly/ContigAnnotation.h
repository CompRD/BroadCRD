// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#ifndef ASSEMBLY_CONTIGANNOTATION_H
#define ASSEMBLY_CONTIGANNOTATION_H

#include "assembly/Annotation.h"
#include "assembly/AssemblyContig.h"
#include "assembly/Location.h"

#include <functional>

class ContigAnnotation
  : public Location<Contig,Annotation>
{
 public:
  ContigAnnotation()
    : Location<Contig,Annotation>()
  {}

  ContigAnnotation( const Location<Contig,Annotation> &other )
    : Location<Contig,Annotation>( other )
  {}

  ContigAnnotation( const Contig &theContig, const Annotation &theAnnotation, 
                    const Interval theInterval )
    : Location<Contig,Annotation>( theContig, theAnnotation, theInterval, orient_FW ) 
  {}

  Contig GetContig( ) const { return m_container; }
  
  const Annotation& GetAnnotation( ) const { return m_element; }

  const String& GetCategory( ) const { return m_element.GetCategory(); }


  // Ordering functors.

  struct orderByContig 
    : public binary_function<ContigAnnotation,ContigAnnotation,bool>
  {
    bool operator() ( const ContigAnnotation& lhs, const ContigAnnotation& rhs ) const
    {
      return lhs.GetContainer() < rhs.GetContainer(); 
    }
  };
};

#endif
