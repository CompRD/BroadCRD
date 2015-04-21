///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// A BasicSequence is a wrapper for a basevector that responds to the
// methods GetBases(), GetQuals() and the output operator.  In its
// constructor, it takes a basevector, an integer id, and a constant
// value that will be used to lazily construct a qualvector for
// GetQuals().

// This is largely intended for use with BoundAlignments, which
// requires the three methods above to be implemented.

#ifndef ASSEMBLY_BASICSEQUENCE_H
#define ASSEMBLY_BASICSEQUENCE_H

#include "Basevector.h"
#include "Qualvector.h"
#include "assembly/BoundAlignment.h"

// TODO: potentially dangerous truncation of index by int id's
class BasicSequence
{
 public:
  BasicSequence( const basevector* p_bases, 
                 const int id = 0,
                 const unsigned char qualScore = 0 )
    : mp_bases( p_bases ),
      mp_quals( 0 ),
      m_id( id ),
      m_qualScore( qualScore )
  {}

  BasicSequence( const BasicSequence& other )
    : mp_bases( other.mp_bases ),
      mp_quals( 0 ),
      m_id( other.m_id ),
      m_qualScore( other.m_qualScore )
  {}

  BasicSequence& operator=( const BasicSequence& other )
  {
    this->mp_bases    = other.mp_bases;
    if ( mp_quals )
      delete mp_quals;
    this->mp_quals    = 0;
    this->m_id        = other.m_id;
    this->m_qualScore = other.m_qualScore;
    
    return *this;
  }

  ~BasicSequence()
  {
    if ( mp_quals )
      delete mp_quals;
  }


  const basevector& GetBases() const { return *mp_bases; }

  const qualvector& GetQuals() const
  {
    if ( mp_quals == 0 )
    {
      mp_quals = new qualvector( mp_bases->size() );
      fill( mp_quals->begin(), mp_quals->end(), m_qualScore );
    }
    return *mp_quals;
  }
  
  friend
  ostream& operator<<( ostream& out, const BasicSequence& seq )
  {
    seq.mp_bases->Print(out, seq.m_id);
    return out;
  }

  int GetId() const { return m_id; }
  
  void SetId( const int id ) { m_id = id; }

 private:
  const basevector* mp_bases;
  mutable qualvector* mp_quals;
  int m_id;
  unsigned char m_qualScore;
};

typedef SequenceIterator<BasicSequence> BasicSequenceIterator;
typedef BoundAlignment<BasicSequence, BasicSequence> BasicBoundAlignment;

#endif
