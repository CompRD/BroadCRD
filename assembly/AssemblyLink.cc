///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "assembly/AssemblyLink.h"

#include "assembly/AssemblyContig.h"
#include "assembly/AssemblyRead.h"

#include <algorithm>

Link::Link()
    : mIsKnowable( false )
{
}

Link::Link( const ContigLocation &contigLoc1, const ReadLocation &readLoc1,
            const ContigLocation &contigLoc2, const ReadLocation &readLoc2,
            const ReadPair &thePair )
    : mContigLoc1( contigLoc1 ), mContigLoc2( contigLoc2 ),
      mReadLoc1( readLoc1 ), mReadLoc2( readLoc2 ), mPair( thePair )
{
    if ( mContigLoc1.GetContig() != mReadLoc1.GetContig() )
    {
        cout << mContigLoc1.GetContig() << " != " << mReadLoc1.GetContig() << endl;
        ForceAssert( mContigLoc1.GetContig() == mReadLoc1.GetContig() );
    }

    if ( mContigLoc2.GetContig() != mReadLoc2.GetContig() )
    {
        cout << mContigLoc2.GetContig() << " != " << mReadLoc2.GetContig() << endl;
        ForceAssert( mContigLoc2.GetContig() == mReadLoc2.GetContig() );
    }

    if ( mReadLoc2 < mReadLoc1 )
    {
        swap( mContigLoc1, mContigLoc2 );
        swap( mReadLoc1, mReadLoc2 );
    }

    // create ReadLocationInSuper objects.
    ReadLocationInSuper readSuperLoc1( mContigLoc1, mReadLoc1 );
    ReadLocationInSuper readSuperLoc2( mContigLoc2, mReadLoc2 );

    // Make sure the read locations refer to reads in the pair.
    ReadToken partner = mPair.GetOtherRead( mReadLoc1.GetRead() );

    ForceAssert( partner == mReadLoc2.GetRead() );

    if ( readSuperLoc1.GetSuper() != readSuperLoc2.GetSuper() )
    {
        mIsKnowable = false;
        mIsLogical = false;
    }
    else
    {
        mIsKnowable = true;
        if ( readSuperLoc1.GetOrientation() == readSuperLoc2.GetOrientation() )
            mIsLogical = false;
        else
        {
            mIsLogical = true;

            // Make sure forward read is in superReadLoc1.
            if ( readSuperLoc1.GetOrientation() == orient_RC )
                swap( readSuperLoc1, readSuperLoc2 );

            mObservedSep = readSuperLoc2.Begin() - readSuperLoc1.End();
            mStretch = mObservedSep - mPair.GetExpectedSep();
            mStretch /= mPair.GetExpectedStdDev();
        }
    }
}

int
Link::GetMeasuredInsertSize() const
{
    Read read1( mReadLoc1.GetRead() );
    Read read2( mReadLoc2.GetRead() );

    return ( read1.GetLeftTrim() + read1.GetLength() +
             this->GetMeasuredSep() +
             read2.GetLength() + read2.GetLeftTrim() );
}

bool
Link::HasOverlapWith( const Link & other ) const
{
     // TODO: only handles within a contig, do for intercontig next
     pair<ReadLocation,ReadLocation> thePair = other.GetReadLocations();
     if ( thePair.first.GetContig() != mReadLoc1.GetContig() )
        {
	cout << thePair.first.GetContig() << " != " << mReadLoc1.GetContig() << endl;
	ForceAssert( thePair.first.GetContig() == mReadLoc1.GetContig() );
	}
     if ( thePair.second.GetContig() != mReadLoc2.GetContig() )
        {
	cout << thePair.second.GetContig() << " != " << mReadLoc2.GetContig() << endl;
	ForceAssert( thePair.second.GetContig() == mReadLoc2.GetContig() );
	}

     Interval otherInterval(thePair.first.Begin(), thePair.second.End());
     Interval thisInterval(mReadLoc1.Begin(), mReadLoc2.End());

     return(thisInterval.HasOverlapWith( otherInterval ));
}


void
BuildLinksFromReadLocation( const ReadLocation &theReadLoc,
                            vec<Link> &links,
                            bool append )
{
  if ( ! append )
    links.clear();

  Read read1 = theReadLoc.GetRead();
  ReadPair thePair = read1.GetPair();

  if ( ! thePair.IsValid() )
    return;

  vec<ReadLocation> read2Locs;
  Read read2 = thePair.GetRead2(); read2.GetLocations( read2Locs );

  vec<ContigLocation> contig1Locs, contig2Locs;

  Contig contig1 = theReadLoc.GetContig();
  contig1.GetSelfLocations( contig1Locs );

  for ( vec<ReadLocation>::iterator read2LocIter = read2Locs.begin();
        read2LocIter != read2Locs.end(); ++read2LocIter )
  {
    Contig contig2 = read2LocIter->GetContig();
    contig2.GetSelfLocations( contig2Locs );

    for ( vec<ContigLocation>::iterator contig1LocIter = contig1Locs.begin();
          contig1LocIter != contig1Locs.end(); ++contig1LocIter )
      for ( vec<ContigLocation>::iterator contig2LocIter = contig2Locs.begin();
            contig2LocIter != contig2Locs.end(); ++contig2LocIter )
      {
        links.push_back( Link( *contig1LocIter, theReadLoc,
                               *contig2LocIter, *read2LocIter,
                               thePair ) );
      }
  }
}


void
BuildLinksFromPair( const ReadPair &thePair,
                    vec<Link> &links,
                    bool append )
{
  if ( ! append )
    links.clear();

  vec<ReadLocation> read1Locs, read2Locs;
  vec<ContigLocation> contig1Locs, contig2Locs;

  Read read1 = thePair.GetRead1(); read1.GetLocations( read1Locs );
  Read read2 = thePair.GetRead2(); read2.GetLocations( read2Locs );

  for ( vec<ReadLocation>::iterator read1LocIter = read1Locs.begin();
        read1LocIter != read1Locs.end(); ++read1LocIter )
  {
    Contig contig1 = read1LocIter->GetContig();
    contig1.GetSelfLocations( contig1Locs );

    for ( vec<ReadLocation>::iterator read2LocIter = read2Locs.begin();
          read2LocIter != read2Locs.end(); ++read2LocIter )
    {
      Contig contig2 = read2LocIter->GetContig();
      contig2.GetSelfLocations( contig2Locs );

      for ( vec<ContigLocation>::iterator contig1LocIter = contig1Locs.begin();
            contig1LocIter != contig1Locs.end(); ++contig1LocIter )
        for ( vec<ContigLocation>::iterator contig2LocIter = contig2Locs.begin();
              contig2LocIter != contig2Locs.end(); ++contig2LocIter )
        {
          links.push_back( Link( *contig1LocIter, *read1LocIter,
                                 *contig2LocIter, *read2LocIter,
                                 thePair ) );
        }
    }
  }
}


