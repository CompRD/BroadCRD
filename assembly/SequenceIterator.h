///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// SequenceIterator is a class that maintains a pointer to a
// "sequence," i.e., something that has the methods GetBases(),
// returning a const basevector reference, and GetQuals(), returning a
// const qualvector reference.  

// A SequenceIterator can be thought of as a const iterator, in that
// it points to a given position in that "sequence," and can return
// the base and qual at that position, but cannot modify that
// information.
//
// A SequenceIterator is a random-access iterator, i.e. it can be
// moved an arbitrary distance forward and backward through the
// sequence in constant time.

#ifndef ASSEMBLY_SEQUENCEITERATOR_H
#define ASSEMBLY_SEQUENCEITERATOR_H


#include "assembly/Orientation.h"

template <class Sequence>
class SequenceIterator
{
 public:
  SequenceIterator( const Sequence &sequence, 
		    Orientation orient = orient_FW, int offset = 0)
    : m_sequence( sequence )
  {
    if ( orient == orient_FW )
    {
      for ( int base = 0; base < 4; ++base )
        m_baseTable[base] = base;
      m_incrementSign = 1;
      m_offset = offset;
    }
    else
    {
      for ( int base = 0; base < 4; ++base )
        m_baseTable[base] = 3 - base;
      m_incrementSign = -1;
      m_offset = m_sequence.GetBases().size() - offset - 1;
    }
  }

  unsigned char GetBase() const;
  unsigned char GetQual() const;

  ///Returns distance from start of iterator.
  ///So for a reverse-oriented iterator it returns distance from the end.
  int GetOffset() const { 
    if (m_incrementSign == 1) return  m_offset;
    else return m_sequence.GetBases().size() - m_offset - 1;
  }

  // does not use GetOffset() to limit calls to GetBases().size().
  bool IsPastEnd() const {
    if (m_incrementSign == -1) return  m_offset < 0;
    else return m_offset >= int(m_sequence.GetBases().size());
  }

  bool IsBeforeStart() const {
    return  GetOffset() < 0;
  }

  const Sequence & GetSequence() const { return m_sequence; }

  SequenceIterator & operator++();      // prefix op
  SequenceIterator   operator++( int ); // postfix op

  SequenceIterator & operator--();      // prefix op
  SequenceIterator   operator--( int ); // postfix op

  SequenceIterator & operator+=( int amount );
  SequenceIterator & operator-=( int amount );

  ///Compare only orientation and offset, not sequences!
  friend bool operator==(const SequenceIterator & lhs,
			 const SequenceIterator & rhs) {
    return ( lhs.m_incrementSign == rhs.m_incrementSign 
	     &&  lhs.m_offset == rhs.m_offset);
  }

  ///Compare only orientation and offset, not sequences!
  friend bool operator!=(const SequenceIterator & lhs,
			 const SequenceIterator & rhs) {
    return !(lhs == rhs);
  }

  ///Return difference between the offsets taking into account orientation.
  /// Useful for example if you want to cut out a basevector between two
  /// iterators: this will give you the size of that basevector.
  /// Note that this works even if the two iterators have different
  /// directions. Not sure if that is a bug or a feature!
  friend int operator-(const SequenceIterator & lhs,
		       const SequenceIterator & rhs) {
    return (lhs.GetOffset() - rhs.GetOffset());
  }

 private:
  void Advance( const int amount );

  Sequence m_sequence;
  int m_offset;

  // Dependent upon orientation.
  char m_baseTable[4];
  int  m_incrementSign;
};
  
template <class Sequence>
inline
SequenceIterator<Sequence> operator+(const SequenceIterator<Sequence> & s, 
				     int i) {
  SequenceIterator<Sequence> ret(s);
  ret += i;
  return ret;
}

template <class Sequence>
inline
SequenceIterator<Sequence> operator-(const SequenceIterator<Sequence> & s, 
				     int i) {
  SequenceIterator<Sequence> ret(s);
  ret -= i;
  return ret;
}

template <class Sequence>
inline
SequenceIterator<Sequence> operator+(int i, const SequenceIterator<Sequence> & s) {
  SequenceIterator<Sequence> ret(s);
  ret += i;
  return ret;
}

template <class Sequence>
inline
SequenceIterator<Sequence> operator-(int i, const SequenceIterator<Sequence> & s) {
  SequenceIterator<Sequence> ret(s);
  ret -= i;
  return ret;
}


// Utility methods to generator SequenceIterators for a given Sequence.

template <class Sequence>
inline
SequenceIterator<Sequence>
GetBeginOfSequence( const Sequence &sequence, Orientation orient )
{
  return SequenceIterator<Sequence>( sequence, orient, 0 );
}


template <class Sequence>
inline
SequenceIterator<Sequence>
GetEndOfSequence( const Sequence &sequence, Orientation orient )
{
  return SequenceIterator<Sequence>( sequence, orient, sequence.GetBases().size() );
}


template <class Sequence>
inline
unsigned char 
SequenceIterator<Sequence>::GetBase() const
{
  return m_baseTable[ m_sequence.GetBases()[ m_offset ] ];
}


template <class Sequence>
inline
unsigned char 
SequenceIterator<Sequence>::GetQual() const
{
  return m_sequence.GetQuals()[ m_offset ];
}

template <class Sequence>
inline
void
SequenceIterator<Sequence>::Advance( const int amount )
{
  m_offset += m_incrementSign * amount;
}

// The operators are implemented here so their calls to Advance() can
// be inlined.

// prefix operator++
template <class Sequence>
inline
SequenceIterator<Sequence> & 
SequenceIterator<Sequence>::operator++()
{
  this->Advance( 1 );
  return *this;
}

// postfix operator++
template <class Sequence>
SequenceIterator<Sequence>
SequenceIterator<Sequence>::operator++( int )
{
  SequenceIterator<Sequence> copy( *this );
  this->Advance( 1 );
  return copy;
}


// prefix operator--
template <class Sequence>
SequenceIterator<Sequence> & 
SequenceIterator<Sequence>::operator--()
{
  this->Advance( -1 );
  return *this;
}

// postfix operator--
template <class Sequence>
SequenceIterator<Sequence>
SequenceIterator<Sequence>::operator--( int )
{
  SequenceIterator<Sequence> copy( *this );
  this->Advance( -1 );
  return copy;
}

// operator +=
template <class Sequence>
SequenceIterator<Sequence> &
SequenceIterator<Sequence>::operator+=( int amount )
{
  this->Advance( amount );
  return *this;
}

// operator -=
template <class Sequence>
SequenceIterator<Sequence> &
SequenceIterator<Sequence>::operator-=( int amount )
{
  this->Advance( -amount );
  return *this;
}

#endif
