///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSEMBLY_INTERVAL
#define ASSEMBLY_INTERVAL

#include <iostream>
#include <algorithm>

using namespace std;

class Interval
{
  public:
    Interval( ) 
        : mBegin(0), mEnd(0) { }
    
    Interval( int begin, int end ) 
        : mBegin(begin), mEnd(end) { }
    
    int Begin( )  const { return mBegin; }
    int End( )    const { return mEnd; }
    int Length( ) const { return mEnd - mBegin; }
    int First( )  const { return mBegin; }
    int Last( )   const { return mEnd - 1; }

    Interval CopyShiftedBy( const int amount ) const
    {
        return Interval( mBegin + amount, mEnd + amount );
    }
    
    Interval CopyExtendedBy( const int amount ) const
    {
        return Interval( mBegin, mEnd + amount );
    }
    
    void SetBegin( const int Begin ) { mBegin = Begin; }
    void SetEnd ( const int End ) { mEnd = End; }
 
    int GetOverlapWith( const Interval & other ) const
    {
	int begin = std::max(mBegin, other.Begin());
	int end = std::min(mEnd, other.End());
	return std::max(0, end - begin);
    }

    bool HasOverlapWith( const Interval & other ) const
    {
	if ( GetOverlapWith( other ) )
	   return true;
	else
	   return false;
    }

    Interval GetIntersectionWith( const Interval &other ) const
    // Returns an interval indicating the subset of this interval that
    // intersects the other interval.  [1000,2000) intersected with
    // [1500,2500) is [1500,2000).  If other begins after this ends,
    // the intersection is [this->End(), this->End()).  If other ends before
    // this begins, the intersection is [this->Begin(), this->Begin()).
    {
      int begin = std::max( mBegin, other.mBegin );
      int end   = std::min( mEnd,   other.mEnd );
      if ( end < begin )
      {
        if ( mEnd <= other.mBegin )
          begin = end = mEnd;
        else
          begin = end = mBegin;
      }
      
      return Interval( begin, end );
    }

    Interval GetCoverageBy( const Interval &other ) const
    // Returns an interval indicating what portion of this interval is
    // intersected by the other interval.  One use of this would be to
    // find which bases of a read are overlapped by another, given
    // their locations in a contig.  This function is not commutative
    // the way GetIntersectionWith() is, e.g., [1000,2000) covered by
    // [1500,2500) is [500,1000), i.e., the last 500 bases of
    // [1000,2000) are covered by [1500,2500), while [1500,2500)
    // covered by [1000,2000) is [0,500), i.e., the first 500 bases of
    // [1500,2500) are covered by [1000,2000).
    {
      return this->GetIntersectionWith( other ).CopyShiftedBy( -mBegin );
    }

    bool Contains( const int point )
    {
	return ( point >= mBegin && point <= mEnd );
    }

    bool operator< ( const Interval &other ) const
    {
        return ( mBegin < other.mBegin ||
                 mBegin == other.mBegin && mEnd < other.mEnd );
    }
    
    bool operator> ( const Interval &other ) const
    {
        return ( mBegin > other.mBegin ||
                 mBegin == other.mBegin && mEnd > other.mEnd );
    }

    bool operator== ( const Interval &other ) const
    {
        return ( mBegin == other.mBegin &&
                 mEnd == other.mEnd );
    }
    
    bool operator!= ( const Interval &other ) const
    {
        return ! ( *this == other );
    }

    friend ostream &operator<< ( ostream &out, const Interval &anInterval )
    {
      out << anInterval.mBegin << ":" << anInterval.mEnd;
      return out;
    }
    
  private:
    int mBegin;
    int mEnd;
};

#endif
