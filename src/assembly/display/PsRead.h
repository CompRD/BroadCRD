// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#ifndef ASSEMBLY_DISPLAY_PSREAD_H
#define ASSEMBLY_DISPLAY_PSREAD_H

#include "Vec.h"

#include "assembly/AssemblyRead.h"
#include "assembly/ReadLocationInContig.h"
#include "assembly/ContigLocation.h"

#include <map>

class PsReadRegistry;

class PsRead
{
 public:
  PsRead( const ReadLocation& theReadLoc,
          const float originX,
          const float originY )
    : m_read( theReadLoc.GetRead() ),
      m_originX( originX ),
      m_originY( originY ),
      m_tailXDiff( ( theReadLoc.IsForward() ? theReadLoc.Length() : -theReadLoc.Length() ) ),
      m_repetitive( Read(theReadLoc.GetRead()).IsRepetitive() )
  {}

  float GetOriginX() const { return m_originX; }
  float GetOriginY() const { return m_originY; }

  bool IsForward() const { return ( m_tailXDiff > 0 ); }
  bool IsReverse() const { return ( ! IsForward() ); }

  void PrintFilled( ostream& out ) const { this->Print( out, true ); }
  void PrintHollow( ostream& out ) const { this->Print( out, false ); }

  void PrintId( ostream &out ) const;
  void PrintName( ostream &out ) const;

  static void PrintHeader( ostream& out );

  static const int s_readHeight = 50;

 private: 
  void Print( ostream &out, const bool filled ) const;

  Read  m_read;
  float m_originX;
  float m_originY;
  float m_tailXDiff;
  Bool  m_repetitive;
};



#endif
