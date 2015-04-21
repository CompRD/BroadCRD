// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#ifndef ASSEMBLY_DISPLAY_PSLINK_H
#define ASSEMBLY_DISPLAY_PSLINK_H

#include "assembly/AssemblyLink.h"

#include "assembly/display/PsRead.h"

// This is not pretty.  If both p_read1 and p_read2 are not null, it
// prints a whole link; otherwise it prints a dashed half-link.  This
// should probably be broken into two classes: PsHalfLink (which can
// be printed whole or dashed) and PsLink (composed of two
// PsHalfLinks).

class PsLink
{
 public:
  PsLink( const PsRead* p_read1,
          const SuperToken& super1,
          const PsRead* p_read2,
          const SuperToken& super2,
          const int expectedSize,
          const int expectedStdDev,
          const float importance  )
    : mp_read1( p_read1 ),
      mp_read2( p_read2 ),
      m_super1( super1 ),
      m_super2( super2 ),
      m_expectedSize( expectedSize ),
      m_expectedStdDev( expectedStdDev ),
      m_importance( importance )
  {
    if ( ! mp_read1 && mp_read2 )
    {
      swap( mp_read1, mp_read2 );
      swap( m_super1, m_super2 );
    }
  }

  void Print( ostream& out ) const;

  static void PrintHeader( ostream& out );

  float GetImportance() const { return m_importance; }

  int GetExpectedSize() const { return m_expectedSize; }

  friend
  bool operator< ( const PsLink& lhs, const PsLink& rhs )
  {
    return ( lhs.m_importance < rhs.m_importance );
  }

 private:
  const PsRead* mp_read1;
  const PsRead* mp_read2;
  SuperToken m_super1;
  SuperToken m_super2;
  int m_expectedSize;
  int m_expectedStdDev;
  float m_importance;
};

#endif
