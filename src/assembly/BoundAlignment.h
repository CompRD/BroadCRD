///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSEMBLY_BOUNDALIGNMENT
#define ASSEMBLY_BOUNDALIGNMENT

#include "math/Arith.h"
#include "PackAlign.h"
#include "PrintAlignment.h"
#include "ScoreAlignment.h"

#include "assembly/Interval.h"
#include "assembly/Orientation.h"
#include "assembly/SequenceIterator.h"

/// A BoundAlignment ties an alignment to the two Sequence objects
/// involved.  It tracks the orientation of both sequences, and can be
/// flipped and reversed.  It can calculate its score and print itself
/// out.  It can also return a (forward) iterator to itself by which it
/// can be traversed.

template <class SequenceType1, class SequenceType2>
class BoundAlignment
{
 public:
  BoundAlignment( ) { }
  
  BoundAlignment( const SequenceType1 &theFirstSequence,
                  const Orientation    theFirstOrientation,
                  const SequenceType2 &theSecondSequence,
                  const Orientation    theSecondOrientation,
                  const align         &theAlign )
    : m_firstSequence( theFirstSequence ),
      m_firstOrient( theFirstOrientation ),
      m_secondSequence( theSecondSequence ),
      m_secondOrient( theSecondOrientation ),
      m_align( theAlign ),
      m_score( -1 ) 
  { }
  
  const SequenceType1 & GetFirstSequence( ) const { return m_firstSequence; }
  const SequenceType2 & GetSecondSequence( ) const { return m_secondSequence; }
  
  Orientation   GetFirstOrientation( )  const { return m_firstOrient; }
  Orientation   GetSecondOrientation( ) const { return m_secondOrient; }
  
  Interval      GetIntervalOnFirst( )   const { return Interval( m_align.pos1(), m_align.Pos1() ); } 
  Interval      GetIntervalOnSecond( )  const { return Interval( m_align.pos2(), m_align.Pos2() ); }
                  
  /// By default, the score is set to -1.  If the score is less than
  /// zero, when GetScore() is called, the score will be computed (and
  /// cached) using ScoreAlignment().
  Float         SetScore( const Float &score );
  Float         GetScore( )             const;

  int           GetOffset()             const { return m_align.pos1() - m_align.pos2(); }

  int           GetLength()             const;  // Length including all gaps.
  int           GetLengthOnFirst( )     const;  // Length spanned on first sequence.
  int           GetLengthOnSecond( )    const;  // Length spanned on second sequence.

 private:
  int           GetLength( bool includeInsertionsOnFirst, 
                           bool includeInsertionsOnSecond ) const;

 public:
  const align & GetAlign( )         const { return m_align; }

  /// Both SequenceType1 and SequenceType2 must have an output operator.
  /// Skip PrintVisualAlignment if brief = true.
  void Print( ostream &out, bool brief = false );

  /// Flip the orientation of the first and second sequences and
  /// reverse the align.
  void ReverseThis( );

  /// Swap sequence1 and sequence2, which may create a BoundAlignment
  /// of a different type.
  BoundAlignment<SequenceType2,SequenceType1> GetFlipped( ) const;

  /// A class that iterates forward through a BoundAlignment.  It
  /// always starts at the beginning of the alignment, and knows
  /// whether it is at the end of the alignment (Iterator::IsAtEnd()).
  ///
  /// BoundAlignment has the convenience method GetIterator().
        
  class Iterator {
    friend class BoundAlignmentTester;
   public:
    Iterator( const BoundAlignment<SequenceType1,SequenceType2> *p_iteratee )
      : mp_iteratee( p_iteratee ),
        m_alignBlockIdx( 0 ),
        m_basesLeftInGap( mp_iteratee->GetAlign().Gaps( 0 ) ),
        m_basesLeftInLen( mp_iteratee->GetAlign().Lengths( 0 ) ),
        m_iterOnFirst( SequenceIterator<SequenceType1>
		       ( mp_iteratee->GetFirstSequence(),
			 mp_iteratee->GetFirstOrientation(),
			 mp_iteratee->GetAlign().pos1() ) 
		       ),
        m_iterOnSecond( SequenceIterator<SequenceType2>
			( mp_iteratee->GetSecondSequence(),
			  mp_iteratee->GetSecondOrientation(),
			  mp_iteratee->GetAlign().pos2() ) 
			)
    { }
    
   private:
    void Advance() 
    {
      Assert(!IsAtEnd());
      if (IsBeforeStart()) { 
	m_alignBlockIdx = 0; 
	//note that basesLeft should be at max for both, and that the
	//sequence iterators should be right at the beginning.
	return; 
      }
	
      if ( m_basesLeftInGap < 0 )
      {
        ++m_iterOnFirst;
        ++m_basesLeftInGap;
      }
      else if ( m_basesLeftInGap > 0 )
      {
        ++m_iterOnSecond;
        --m_basesLeftInGap;
      }
      else if ( m_basesLeftInLen > 0 )
      {
        ++m_iterOnFirst;
        ++m_iterOnSecond;
        --m_basesLeftInLen;
      }
      
      if ( m_basesLeftInGap == 0 &&
           m_basesLeftInLen == 0 )
        if ( ++m_alignBlockIdx < mp_iteratee->GetAlign().Nblocks() )
        {
          m_basesLeftInGap = mp_iteratee->GetAlign().Gaps( m_alignBlockIdx );
          m_basesLeftInLen = mp_iteratee->GetAlign().Lengths( m_alignBlockIdx );
        }
        else
          m_alignBlockIdx = mp_iteratee->GetAlign().Nblocks();
    }

    void GoBack()  {
      Assert(!IsBeforeStart());
      if (IsAtEnd()) { 
	--m_alignBlockIdx;
	m_basesLeftInGap = 0;
	m_basesLeftInLen = 1;
	--m_iterOnFirst;
	--m_iterOnSecond;
	return;
      }

      int len = mp_iteratee->GetAlign().Lengths( m_alignBlockIdx);
      int gap = mp_iteratee->GetAlign().Gaps( m_alignBlockIdx );
      //Now update our pointers if we got to the beginning of a block.
      if ( m_basesLeftInGap == gap 
	   && m_basesLeftInLen == len ) {
	if ( --m_alignBlockIdx >=0 ) {
	  m_basesLeftInGap = 0;
	  m_basesLeftInLen = 0;
	}
      }
      // First, decrease the position if so indicated.
      if ( m_basesLeftInLen !=  len) {
	--m_iterOnFirst;
	--m_iterOnSecond;
	++m_basesLeftInLen;
      } 
      else if (  gap < 0  && m_basesLeftInGap != gap) {
	--m_iterOnFirst;
	--m_basesLeftInGap;
      }
      else if ( gap > 0 && m_basesLeftInGap != gap ) {
	--m_iterOnSecond;
	++m_basesLeftInGap;
      }
    }

   public:
    const BoundAlignment * Alignment() const { return mp_iteratee; }

    bool IsAtEnd() { 
      return m_alignBlockIdx == mp_iteratee->GetAlign().Nblocks(); 
    }
    bool IsBeforeStart() { return m_alignBlockIdx < 0; }
      

    Iterator & operator++()
    {
      this->Advance();
      return *this;
    }

    Iterator operator++(int)
    {
      Iterator copy( *this );
      this->Advance();
      return copy;
    }

    Iterator & operator--() {
      this->GoBack();
      return *this;
    }

    Iterator operator--(int) {
      Iterator copy( *this );
      this->GoBack();
      return copy;
    }

    bool IsGap() const         { return ( m_basesLeftInGap != 0 ); }

    ///Return true if we are between bases on the first read.
    bool IsGapOnFirst() const  { return ( m_basesLeftInGap > 0 ); }

    ///Return true if we are between bases on the second read.
    bool IsGapOnSecond() const { return ( m_basesLeftInGap < 0 ); }
    
    bool IsMatch() const { 
      return ( ! IsGap() 
	       && m_iterOnFirst.GetBase() == m_iterOnSecond.GetBase() ); 
    }

    ///Returns length of gap (+ve for gap on 2nd, -ve for gap on 1st), or false if not in gap.
    int GetGapLength() const { return IsGap() ? mp_iteratee->GetAlign().Gaps( m_alignBlockIdx ) : false; }

    int           GetOffsetOnFirst() const { return m_iterOnFirst.GetOffset(); }
    unsigned char GetBaseOnFirst()   const { return m_iterOnFirst.GetBase(); }
    unsigned char GetQualOnFirst()   const { return m_iterOnFirst.GetQual(); }

    const SequenceIterator<SequenceType1> &  GetIterOnFirst()  { return m_iterOnFirst; }
    
    int           GetOffsetOnSecond() const { return m_iterOnSecond.GetOffset(); }
    unsigned char GetBaseOnSecond()   const { return m_iterOnSecond.GetBase(); }
    unsigned char GetQualOnSecond()   const { return m_iterOnSecond.GetQual(); }
    
    const SequenceIterator<SequenceType2> &  GetIterOnSecond() { return m_iterOnSecond; }

    const SequenceType1 & GetFirst() const {
      return mp_iteratee->GetFirstSequence();
    }
    const SequenceType2 & GetSecond() const {
      return mp_iteratee->GetSecondSequence();
    }

    ///Does not check to see if the underlying BoundAlignments are the same.
    friend bool operator==(const Iterator & l, const Iterator & r) {
      return l.m_alignBlockIdx == r.m_alignBlockIdx
	&& l.m_basesLeftInGap ==r.m_basesLeftInGap
	&& l.m_basesLeftInLen == r.m_basesLeftInLen 
	&& l.GetOffsetOnFirst() == r.GetOffsetOnFirst()
	&& l.GetOffsetOnSecond() == r.GetOffsetOnSecond();
    }

    ///Does not check to see if the underlying BoundAlignments are the same.
    friend bool operator!=(const Iterator & l, const Iterator & r) {
      return !(l==r);
    }

  private:
    const BoundAlignment<SequenceType1,SequenceType2> *mp_iteratee;
    int m_alignBlockIdx;
    int m_basesLeftInGap;
    int m_basesLeftInLen;
    SequenceIterator<SequenceType1> m_iterOnFirst;
    SequenceIterator<SequenceType2> m_iterOnSecond;
  };

  Iterator GetIterator() const { return Iterator( this ); }

 private:
  SequenceType1 m_firstSequence;
  Orientation   m_firstOrient;
  SequenceType2 m_secondSequence;
  Orientation   m_secondOrient;
  align         m_align;
  mutable Float m_score;
};
  


// Implementations


template <class SequenceType1, class SequenceType2>
inline
int
BoundAlignment<SequenceType1,SequenceType2>::GetLength( bool includeInsertionsOnFirst,
                                                        bool includeInsertionsOnSecond ) const
{
  int sumLengths = 0;
  int sumGaps = 0;
  
  for ( int blockIdx = 0; blockIdx < m_align.Nblocks(); ++blockIdx )
  {
    int gap = m_align.Gaps( blockIdx );
    int len = m_align.Lengths( blockIdx );
    if ( gap < 0 && includeInsertionsOnFirst )
      sumGaps -= gap;
    if ( gap > 0 && includeInsertionsOnSecond )
      sumGaps += gap;
    sumLengths += len;
  }

  return sumLengths + sumGaps;
}


template <class SequenceType1, class SequenceType2>
inline
int
BoundAlignment<SequenceType1,SequenceType2>::GetLength( ) const
{
  return this->GetLength( true, true );
}


template <class SequenceType1, class SequenceType2>
inline
int
BoundAlignment<SequenceType1,SequenceType2>::GetLengthOnFirst( ) const
{
  return this->GetLength( true, false );
}


template <class SequenceType1, class SequenceType2>
inline
int
BoundAlignment<SequenceType1,SequenceType2>::GetLengthOnSecond( ) const
{
  return this->GetLength( false, true );
}


template <class SequenceType1, class SequenceType2>
inline
Float
BoundAlignment<SequenceType1,SequenceType2>::GetScore( ) const
{
  if ( m_score < Float(0.0f) )
  {
    static basevector bases1, bases2;
    bases1 = m_firstSequence.GetBases();
    bases2 = m_secondSequence.GetBases();
    static qualvector quals1, quals2;
    quals1 = m_firstSequence.GetQuals();
    quals2 = m_secondSequence.GetQuals();
    
    if ( m_firstOrient == orient_RC )
    {
      bases1.ReverseComplement();
      quals1.ReverseMe();
    }

    if ( m_secondOrient == orient_RC )
    {
      bases2.ReverseComplement();
      quals2.ReverseMe();
    }

    m_score = ScoreAlignment( m_align, bases1, quals1, bases2, quals2 );
  }
  
  return m_score;
}


template <class SequenceType1, class SequenceType2>
inline
void
BoundAlignment<SequenceType1,SequenceType2>::Print( ostream &out, bool brief )
{
    static basevector bases1, bases2;
    bases1 = m_firstSequence.GetBases();
    bases2 = m_secondSequence.GetBases();
    static qualvector quals1, quals2;
    quals1 = m_firstSequence.GetQuals();
    quals2 = m_secondSequence.GetQuals();
    
    if ( m_firstOrient == orient_RC )
    {
      bases1.ReverseComplement();
      quals1.ReverseMe();
    }

    if ( m_secondOrient == orient_RC )
    {
      bases2.ReverseComplement();
      quals2.ReverseMe();
    }

    out << "Alignment between " 
        << m_firstSequence << ":" << m_align.pos1() << "-" << m_align.Pos1() 
	<< "_" << bases1.size( )
        << "(" << m_firstOrient << ")"
        << " and " 
        << m_secondSequence << ":" << m_align.pos2() << "-" << m_align.Pos2() 
	<< "_" << bases2.size( )
        << "(" << m_secondOrient << ")" 
        << endl;

    if ( brief )
      return;

    PrintVisualAlignment( True, out, bases1, bases2, m_align, quals1, quals2 );
}


template <class SequenceType1, class SequenceType2>
inline
void
BoundAlignment<SequenceType1,SequenceType2>::ReverseThis( )
{
    m_align.ReverseThis( m_firstSequence.GetLength(), m_secondSequence.GetLength() );
    m_firstOrient = Flip( m_firstOrient );
    m_secondOrient = Flip( m_secondOrient );
}

template <class SequenceType1, class SequenceType2>
inline
BoundAlignment<SequenceType2,SequenceType1>
BoundAlignment<SequenceType1,SequenceType2>::GetFlipped( ) const
{
    align flippedAlign( m_align );
    flippedAlign.Flip();
    
    return BoundAlignment<SequenceType2,SequenceType1>( m_secondSequence,
                                                        m_secondOrient,
                                                        m_firstSequence,
                                                        m_firstOrient,
                                                        flippedAlign );
}

#endif

