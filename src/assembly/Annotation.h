// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#ifndef ASSEMBLY_ANNOTATION_H
#define ASSEMBLY_ANNOTATION_H

#include "String.h"
#include "Vec.h"

#include <map>

class AnnotationRegistry
{
 private:
  AnnotationRegistry() {}
  
  static AnnotationRegistry* sp_singleton;
  
 public:
  static AnnotationRegistry* GetRegistryPtr()
  {
    if ( sp_singleton == 0 )
    {
      sp_singleton = new AnnotationRegistry();
      // Category "unknown" always has code 0.
      sp_singleton->GetCategoryCode( "unknown" );
    }
    return sp_singleton;
  }
  
  int GetCategoryCode( const String& category )
  {
    map<String,int>::const_iterator lookupIter = m_codeLookup.find( category );
    if ( lookupIter != m_codeLookup.end() &&
	 lookupIter->first == category )
      return lookupIter->second;
    else
    {
      int newCode = m_codeNames.size();
      m_codeNames.push_back( category );
      m_codeLookup.insert( make_pair( category, newCode ) );
      return newCode;
    }
  }
  
  const String& GetCategoryName( const int code )
  {
    if ( code > 0 && code < (int)m_codeNames.size() )
      return m_codeNames[ code ];
    else
      return m_codeNames[0];
  }
  
 private:
  vec<String> m_codeNames;
  map<String,int> m_codeLookup;
};
      

class Annotation
{
 public:
  Annotation()
    : m_code(-1)
  {}

  Annotation( const String& category )
  {
    m_code = AnnotationRegistry::GetRegistryPtr()->GetCategoryCode( category );
  }

  const String& GetCategory() const
  {
    return AnnotationRegistry::GetRegistryPtr()->GetCategoryName( m_code );
  }

  bool operator< ( const Annotation& other ) const
  {
    return m_code < other.m_code;
  }

 private:
  int m_code;
};
  

#endif
