///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// An assessment of a ReadPair where both Reads in the ReadPair are placed.

#ifndef ASSEMBLY_LINK
#define ASSEMBLY_LINK

#include <utility>

#include "assembly/ContigLocation.h"
#include "assembly/ReadLocationInContig.h"
#include "assembly/ReadLocationInSuper.h"
#include "assembly/ReadPair.h"

#include "Vec.h"

#include <math.h>

class Link
{
  public:
    Link();

    Link( const ContigLocation &contigLoc1, const ReadLocation &readLoc1, 
          const ContigLocation &contigLoc2, const ReadLocation &readLoc2, 
          const ReadPair &thePair);

    // These first and second elements of these pairs all correspond,
    // i.e. the first element of each refers to one read, and the
    // second element of each refers to the other.

    // Upone construction, the order of the elements is determined by
    // the ReadLocations.  The ReadLocation of the first read is no
    // greater than that of the second.

    pair<ReadToken,ReadToken>                     GetReads() const;
    pair<ContigToken,ContigToken>                 GetContigs() const;
    pair<SuperToken,SuperToken>                   GetSupers() const;

    pair<ReadLocation,ReadLocation>               GetReadLocations() const;
    pair<ContigLocation,ContigLocation>           GetContigLocations() const;
    pair<ReadLocationInSuper,ReadLocationInSuper> GetReadLocationsInSupers() const;

    ReadToken            GetRead1() const { return mReadLoc1.GetRead(); }
    ReadToken            GetRead2() const { return mReadLoc2.GetRead(); }

    ContigToken          GetContig1() const { return mReadLoc1.GetContig(); }
    ContigToken          GetContig2() const { return mReadLoc2.GetContig(); }

    SuperToken           GetSuper1() const { return mContigLoc1.GetSuper(); }
    SuperToken           GetSuper2() const { return mContigLoc2.GetSuper(); }

    ReadLocation         GetReadLocation1() const { return mReadLoc1; }
    ReadLocation         GetReadLocation2() const { return mReadLoc2; }

    ReadLocationInSuper  GetReadLocationInSuper1() const;
    ReadLocationInSuper  GetReadLocationInSuper2() const;

    ReadPair GetPair() const { return mPair; }

    // Swap the two reads internally.
    void SwapReads() 
    {
      swap( mContigLoc1, mContigLoc2 );
      swap( mReadLoc1, mReadLoc2 );
    }

    // Get the separation observed in the assembly of the read pair.
    // Will assert if the link is not logical.
    int GetMeasuredSep() const; 

    // Get the insert size observed in the assembly of the read pair,
    // inclusive of read lengths and trims.
    // Will assert if the link is not logical.
    int GetMeasuredInsertSize() const; 

    // Get expected separation of the read pair.
    int GetExpectedSep() const;

    // Get expected insert size of the read pair, inclusive of read
    // lengths and trims.
    int GetExpectedInsertSize() const;

    // Get the difference between measured and expected separation.
    // Will assert if link is not logical.
    int GetDiffSep() const; 
  
    // Get the expected standard deviation of the pair's insert size
    // (and therefore separation).
    int GetExpectedStdDev()  const;

    // True if and only if the reads are in different supers.
    bool IsUnknown() const;

    // True if and only if the reads are in same super and are oriented in
    // different directions.
    bool IsLogical() const;

    // If the link is logical, return the result of the calculation 
    //
    //    ( measured separation - expected separation )
    //    ---------------------------------------------
    //            expected standard deviation
    //
    // If the link is not logical, assert.
    float GetStretch() const;

    // True if and only if the link is logical and the absolute value
    // of its stretch is less than or equal to the given threshold.
    bool IsValid( const float threshold = 3.0 ) const;

    // True if and only if the link is logical and its stretch is
    // greater than or equal to lowThreshold and less than or equal to
    // highThreshold.
    bool IsValid( const float lowThreshold, const float highThreshold ) const;

    // True if and only if the other link overlaps this link. Asserts if they are on
    // different contigs.
    bool HasOverlapWith( const Link & other ) const;

    bool InvolvesContig( ContigToken theContig ) const
    {  return ( mContigLoc1.GetContig() == theContig ||
                mContigLoc2.GetContig() == theContig );  }

    bool InvolvesSuper( SuperToken theSuper ) const
    {  return ( mContigLoc1.GetSuper() == theSuper ||
                mContigLoc2.GetSuper() == theSuper );  }

    bool IsBetweenContigs() const
    {  return ! (mContigLoc1 == mContigLoc2);  }

    bool IsBetweenSupers() const
    { 
        return ( mContigLoc1.GetSuper().IsValid() &&
                 mContigLoc2.GetSuper().IsValid() &&
                 mContigLoc1.GetSuper() != mContigLoc2.GetSuper() );
    }

    bool operator== ( const Link &other ) const
    {
        return ( mContigLoc1 == other.mContigLoc1 &&
                 mContigLoc2 == other.mContigLoc2 &&
                 mReadLoc1   == other.mReadLoc1 &&
                 mReadLoc2   == other.mReadLoc2 );
    }

    friend
    ostream & operator<< ( ostream &out, const Link &aLink )
    {
        return out << aLink.mReadLoc1 << "\t" << aLink.mContigLoc1
                   << "\t<=(" << aLink.mPair.GetExpectedSep()
                   << "~" << aLink.mPair.GetExpectedStdDev()
                   << ")=>\t"
                   << aLink.mReadLoc2 << "\t" << aLink.mContigLoc2;
    }

  private:
    ContigLocation mContigLoc1, mContigLoc2;
    ReadLocation mReadLoc1, mReadLoc2;
    ReadPair mPair;
    bool mIsKnowable;
    bool mIsLogical;
    int mObservedSep;
    float mStretch;
};


// BuildLinksFromPair is a utility function that fills out the "links"
// parameter with all the links resulting from all the placements of
// the reads in "thePair".

void
BuildLinksFromPair( const ReadPair &thePair,
                    vec<Link> &links,
                    bool append = false );


// BuildLinksFromReadLocation is a utility function that fills out the
// "links" parameter with all the links resulting from the given read
// placement and all the placements of its partner.

void
BuildLinksFromReadLocation( const ReadLocation &theReadLoc,
                            vec<Link> &links,
                            bool append = false );


// Inlined methods.

inline
pair<ReadToken,ReadToken> 
Link::GetReads() const
{
    return make_pair( mReadLoc1.GetRead(), mReadLoc2.GetRead() );
}

inline
pair<ContigToken,ContigToken> 
Link::GetContigs() const
{
    return make_pair( mContigLoc1.GetContig(), mContigLoc2.GetContig() );
}

inline
pair<SuperToken,SuperToken> 
Link::GetSupers() const
{
    return make_pair( mContigLoc1.GetSuper(), mContigLoc2.GetSuper() );
}

inline
pair<ReadLocation,ReadLocation> 
Link::GetReadLocations() const
{
    return make_pair( mReadLoc1, mReadLoc2 );
}

inline
pair<ContigLocation,ContigLocation> 
Link::GetContigLocations() const
{
    return make_pair( mContigLoc1, mContigLoc2 );
}

inline
ReadLocationInSuper 
Link::GetReadLocationInSuper1() const
{
    return ReadLocationInSuper(mContigLoc1,mReadLoc1);
}

inline
ReadLocationInSuper 
Link::GetReadLocationInSuper2() const
{
    return ReadLocationInSuper(mContigLoc2,mReadLoc2);
}

inline
pair<ReadLocationInSuper,ReadLocationInSuper> 
Link::GetReadLocationsInSupers() const
{
    return make_pair( this->GetReadLocationInSuper1(),
                      this->GetReadLocationInSuper2() );
}

inline
int
Link::GetMeasuredSep() const
{
    ForceAssert( this->IsLogical() );

    return mObservedSep;
}

inline
int
Link::GetExpectedSep() const
{
    return mPair.GetExpectedSep();
}

inline
int
Link::GetExpectedInsertSize() const
{
    return mPair.GetExpectedInsertSize();
}

inline
int
Link::GetDiffSep() const
{
    ForceAssert( this->IsLogical() );

    return (this->GetMeasuredSep() - this->GetExpectedSep());
}

inline
int
Link::GetExpectedStdDev() const
{
    return mPair.GetExpectedStdDev();
}

inline
bool
Link::IsUnknown() const
{
    return ! mIsKnowable;
}

inline
bool
Link::IsLogical() const
{
    return mIsLogical;
}

inline
float 
Link::GetStretch() const
{
    ForceAssert( this->IsLogical() );

    return mStretch;
}

inline
bool 
Link::IsValid( const float threshold ) const
{
    return ( this->IsLogical() && fabs( mStretch ) <= threshold );
}

inline
bool
Link::IsValid( const float lowThreshold, const float highThreshold ) const
{
    float absStretch = fabs( mStretch );
    return ( this->IsLogical() && 
             lowThreshold <= absStretch &&
             absStretch <= highThreshold );
}

#endif
