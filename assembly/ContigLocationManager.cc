///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "assembly/ContigLocationManager.h"

#include "assembly/ContigDataManager.h"
#include "assembly/SuperDataManager.h"

ContigLocationManager::ContigLocationManager( SuperDataManager * &rpSuperDataMgr,
                                              ContigDataManager * &rpContigDataMgr )
    : mrpSuperDataMgr( rpSuperDataMgr ),
      mrpContigDataMgr( rpContigDataMgr ),
      mbDirty( true ),
      mbModified( false ),
      mNumEmpty( 0 )
{
}


void
ContigLocationManager::Add( const ContigLocation &theContigLoc )
{
    if ( mbDirty ) this->Clean();

    int locIdx = mVecContigLocations.size();

    mVecContigLocations.push_back( theContigLoc );

    int contigId = theContigLoc.GetContig().GetId();
    ForceAssert( theContigLoc.GetContig().IsValid() );
    if ( contigId >= static_cast<int>( mIndicesByContig.size() ) )
    {
      if ( contigId >= static_cast<int>( mIndicesByContig.capacity() ) )
        mIndicesByContig.reserve( (contigId + 1) * 2 );
      mIndicesByContig.resize( contigId + 1 );
    }
    mIndicesByContig[ contigId ].push_back( locIdx );

    int superId = theContigLoc.GetSuper().GetId();
    ForceAssert( theContigLoc.GetSuper().IsValid() );
    if ( superId >= static_cast<int>( mIndicesBySuper.size() ) )
    {
      if ( superId >= static_cast<int>( mIndicesBySuper.capacity() ) )
        mIndicesBySuper.reserve( (superId + 1) * 2 );
      mIndicesBySuper.resize( superId + 1 );
    }
    mIndicesBySuper[ superId ].push_back( locIdx );

    if ( superId >= static_cast<int>( mVecSuperBegin.size() ) )
    {
      if ( superId >= static_cast<int>( mVecSuperBegin.capacity() ) )
        mVecSuperBegin.reserve( (superId + 1) * 2 );
      mVecSuperBegin.resize( superId + 1, 0 );
    }

    if ( mIndicesBySuper[ superId ].size() == 1 )
      mVecSuperBegin[ superId ] = theContigLoc.Begin();
    else
      mVecSuperBegin[ superId ] = min( mVecSuperBegin[ superId ], theContigLoc.Begin() );

    if ( superId >= static_cast<int>( mVecSuperEnd.size() ) )
    {
      if ( superId >= static_cast<int>( mVecSuperEnd.capacity() ) )
        mVecSuperEnd.reserve( (superId + 1) * 2 );
      mVecSuperEnd.resize( superId + 1, 0 );
    }

    if ( mIndicesBySuper[ superId ].size() == 1 )
      mVecSuperEnd[ superId ] = theContigLoc.End();
    else
      mVecSuperEnd[ superId ] = max( mVecSuperEnd[ superId ], theContigLoc.End() );

    if ( superId >= static_cast<int>( mSumContigLengths.size() ) ) {
      if ( superId >= static_cast<int>( mSumContigLengths.capacity() ) )
        mSumContigLengths.reserve( (superId + 1) * 2 );
      mSumContigLengths.resize( superId + 1, 0 );
    }
    
    mSumContigLengths[superId] += theContigLoc.Length();

    mbModified = true;
}

bool
ContigLocationManager::Remove( const ContigLocation &theContigLoc )
{
    if ( mbDirty ) Clean();
    
    int locIdx = -1;

    Contig theContig = theContigLoc.GetContig();

    vec<int> &vecIndicesFromContig = mIndicesByContig[ theContig.GetId() ];

    vec<int>::iterator idxIter;
    idxIter = vecIndicesFromContig.begin(); 
    while ( idxIter != vecIndicesFromContig.end() )
    {
        ContigLocation &theLoc = mVecContigLocations[ *idxIter ];

        if ( theContigLoc == theLoc )
        {
            locIdx = *idxIter;
            vecIndicesFromContig.erase( idxIter );
            mbModified = true;
            break;
        }
        else
            ++idxIter;
    }

    if ( locIdx < 0 )
        return false;


    Super theSuper = theContigLoc.GetSuper();

    vec<int> &vecIndicesFromSuper = mIndicesBySuper[ theSuper.GetId() ];

    bool bMustRecalcBeginEnd = false;

    idxIter = vecIndicesFromSuper.begin(); 
    while ( idxIter != vecIndicesFromSuper.end() )
    {
        if ( *idxIter == locIdx )
        {
            vecIndicesFromSuper.erase( idxIter );
            mbModified = true;
            if ( mVecSuperBegin[ theSuper.GetId() ] >= theContigLoc.Begin() ||
                 mVecSuperEnd  [ theSuper.GetId() ] <= theContigLoc.End() )
              bMustRecalcBeginEnd = true;
            break;
        }
        else
          ++idxIter;
    }

    mSumContigLengths[ theSuper.GetId() ] -= theContigLoc.Length();
    
    if ( bMustRecalcBeginEnd )
    {
      int begin = 0;
      int end   = 0;

      for ( idxIter = vecIndicesFromSuper.begin();
            idxIter != vecIndicesFromSuper.end();
            ++idxIter )
      {
        ContigLocation &theLoc = mVecContigLocations[ *idxIter ];

        begin = min<int>( begin, theLoc.Begin() );
        end   = max<int>( end,   theLoc.End() );
      }

      mVecSuperBegin[ theSuper.GetId() ] = begin;
      mVecSuperEnd  [ theSuper.GetId() ] = end;
    }

    mVecContigLocations[ locIdx ] = ContigLocation();
    ++mNumEmpty;

    unsigned int numContigs = mIndicesByContig.size();

    if ( mNumEmpty > numContigs ) {
      // cout << "Purging empty ContigLocations." << endl;
      unsigned int keep = 0;
      for ( unsigned int i = 0; i < mVecContigLocations.size(); ++i )
        if ( mVecContigLocations[i].GetSuper().IsValid() ) {
          if ( keep != i )
            mVecContigLocations[keep] = mVecContigLocations[i];
          ++keep;
        }
      mVecContigLocations.resize( keep );
      mNumEmpty = 0;
      mbDirty = true;
    }

    return true;
}


void
ContigLocationManager::ClearSuper( const int superId )
{
    if ( mbDirty ) Clean();

    // If the super is empty and has been added since we last
    // Clean()ed, the index vector will be too small.  Check that the
    // id is legal and return.
    if ( superId >= (int) mIndicesBySuper.size() )
    {
      ForceAssertLt( superId, mrpSuperDataMgr->GetSize() );
      return;
    }

    // We could use Remove() to do this, but it's rather inefficient for large
    // supers, so we'll do it "by hand".
    for ( vec<int>::iterator superIdxIter = mIndicesBySuper[ superId ].begin();
          superIdxIter != mIndicesBySuper[ superId ].end(); ++superIdxIter ) {

      ContigLocation theContigLoc = mVecContigLocations[ *superIdxIter ];

      Contig theContig = theContigLoc.GetContig();

      vec<int> &vecIndicesFromContig = mIndicesByContig[ theContig.GetId() ];

      for ( vec<int>::iterator contigIdxIter = vecIndicesFromContig.begin(); 
            contigIdxIter != vecIndicesFromContig.end(); ++contigIdxIter ) {
        if ( *contigIdxIter == *superIdxIter ) {
          vecIndicesFromContig.erase( contigIdxIter );
          mbModified = true;
          break;
        }
      }

      mVecContigLocations[ *superIdxIter ] = ContigLocation();
      ++mNumEmpty;
    }

    mIndicesBySuper[ superId ].clear();

    mVecSuperBegin[ superId ] = 0;
    mVecSuperEnd  [ superId ] = 0;
    mSumContigLengths[ superId ] = 0;

    unsigned int numContigs = mIndicesByContig.size();

    if ( mNumEmpty > numContigs ) {
      double startTime = WallClockTime();
      unsigned int keep = 0;
      for ( unsigned int i = 0; i < mVecContigLocations.size(); ++i )
        if ( mVecContigLocations[i].GetSuper().IsValid() ) {
          if ( keep != i )
            mVecContigLocations[keep] = mVecContigLocations[i];
          ++keep;
        }
      mVecContigLocations.resize( keep );
      mNumEmpty = 0;
      mbDirty = true;

      // cout << "Purged empty ContigLocations in " << TimeSince( startTime ) << "." << endl;
    }
}


void
ContigLocationManager::Reverse( const SuperToken theSuperToken )
{
    if ( mbDirty ) Clean();

    // If the super is empty and has been added since we last
    // Clean()ed, the index vector will be too small.  Check that the
    // id is legal and return.
    const int superId = theSuperToken.GetId();

    if ( superId >= (int) mIndicesBySuper.size() )
    {
      ForceAssertLt( superId, mrpSuperDataMgr->GetSize() );
      return;
    }

    vec<int> &vecIndices = mIndicesBySuper[ superId ];

    Super theSuper(theSuperToken);

    int length = theSuper.GetLength();

    for ( vec<int>::iterator idxIter = vecIndices.begin();
          idxIter != vecIndices.end(); ++idxIter )
    {
        ContigLocation &theLoc = mVecContigLocations[ *idxIter ];

        theLoc = ContigLocation( theSuper, theLoc.GetContig(),
                                 Interval( length - theLoc.End(), length - theLoc.Begin() ),
                                 Flip( theLoc.GetOrientation() ) );
    }
}

    
void
ContigLocationManager::Shift( const int superId, const int shiftAmount )
{
    if ( mbDirty ) Clean();

    // If the super is empty and has been added since we last
    // Clean()ed, the index vector will be too small.  Check that the
    // id is legal and return.
    if ( superId >= (int) mIndicesBySuper.size() )
    {
      ForceAssertLt( superId, mrpSuperDataMgr->GetSize() );
      return;
    }

    vec<int> &vecIndices = mIndicesBySuper[ superId ];

    mVecSuperBegin[ superId ] += shiftAmount;
    mVecSuperEnd  [ superId ] += shiftAmount;

    for ( vec<int>::iterator idxIter = vecIndices.begin(); 
          idxIter != vecIndices.end(); ++idxIter )
    {
        ContigLocation &theLoc = mVecContigLocations[ *idxIter ];
        theLoc = ContigLocation( theLoc.GetSuper(), theLoc.GetContig(),
                                 theLoc.GetInterval().CopyShiftedBy( shiftAmount ),
                                 theLoc.GetOrientation() );
    }
}


void
ContigLocationManager::Normalize( const int superId )
{
    if ( mbDirty ) Clean();

    // If the super is empty and has been added since we last
    // Clean()ed, the index vector will be too small.  Check that the
    // id is legal and return.
    if ( superId >= (int) mIndicesBySuper.size() )
    {
      ForceAssertLt( superId, mrpSuperDataMgr->GetSize() );
      return;
    }

    this->Shift( superId, -1 * mVecSuperBegin[ superId ] );
}


int
ContigLocationManager::GetNumContigLocationsInSuper( const int superId )
{
    if ( mbDirty ) Clean();

    const vec<int> &vecIndices = mIndicesBySuper[ superId ];
    return vecIndices.size();
}

void
ContigLocationManager::GetBySuper( const int superId,
                                   vec<ContigLocation> &vecContigLocs )
{
    if ( mbDirty ) Clean();

    vecContigLocs.clear();

    // If the super is empty and has been added since we last
    // Clean()ed, the index vector will be too small.  Check that the
    // id is legal and return.
    if ( superId >= (int) mIndicesBySuper.size() )
    {
      ForceAssertLt( superId, mrpSuperDataMgr->GetSize() );
      return;
    }

    vec<int> &vecIndices = mIndicesBySuper[ superId ];

    vecContigLocs.reserve( vecIndices.size() );

    for ( vec<int>::iterator idxIter = vecIndices.begin(); 
          idxIter != vecIndices.end(); ++idxIter )
        vecContigLocs.push_back( mVecContigLocations[ *idxIter ] );
}


void
ContigLocationManager::GetByContig( const int contigId, 
                                    vec<ContigLocation> &vecContigLocs )
{
    if ( mbDirty ) Clean();

    vecContigLocs.clear();

    // If the contig is not placed and has been added since we last
    // Clean()ed, the index vector will be too small.  Check that the
    // id is legal and return.
    if ( contigId >= (int) mIndicesByContig.size() )
    {
      ForceAssertLt( contigId, mrpContigDataMgr->GetSize() );
      return;
    }

    vec<int> &vecIndices = mIndicesByContig[ contigId ];

    vecContigLocs.reserve( vecIndices.size() );

    for ( vec<int>::iterator idxIter = vecIndices.begin(); 
          idxIter != vecIndices.end(); ++idxIter )
        vecContigLocs.push_back( mVecContigLocations[ *idxIter ] );
}

int
ContigLocationManager::GetSuperLength( const int superId )
{
    if ( mbDirty ) Clean();

    // If the super is empty and has been added since we last
    // Clean()ed, the cached begin vector will be too small.  Check
    // that the id is legal and return zero.
    if ( superId >= (int) mVecSuperBegin.size() )
    {
      ForceAssertLt( superId, mrpSuperDataMgr->GetSize() );
      return 0;
    }

    return ( mVecSuperEnd[ superId ] - mVecSuperBegin[ superId ] );
}

longlong
ContigLocationManager::GetSumOfContigLengths( const int superId )
{
    if ( mbDirty ) Clean();

    if ( superId >= (int) mSumContigLengths.size() )
    {
      ForceAssertLt( superId, mrpSuperDataMgr->GetSize() );
      return 0;
    }

    return mSumContigLengths[ superId ];
}

void 
ContigLocationManager::Clean()
{
    // cout << "Indexing ContigLocationManager." << endl;

    double startTime = WallClockTime();

    mIndicesBySuper.resize( mrpSuperDataMgr->GetSize() );

    for ( unsigned int i = 0; i < mIndicesBySuper.size(); ++i )
      mIndicesBySuper[i].clear();

    for ( unsigned int locIdx = 0; locIdx < mVecContigLocations.size(); ++locIdx )
    {
        int superId = mVecContigLocations[ locIdx ].GetSuper().GetId();

        mIndicesBySuper[ superId ].push_back( locIdx );
    }


    mIndicesByContig.resize( mrpContigDataMgr->GetSize() );

    for ( unsigned int i = 0; i < mIndicesByContig.size(); ++i )
      mIndicesByContig[i].clear();

    for ( unsigned int locIdx = 0; locIdx < mVecContigLocations.size(); ++locIdx )
    {
        int contigId = mVecContigLocations[ locIdx ].GetContig().GetId();

        mIndicesByContig[ contigId ].push_back( locIdx );
    }


    // Cache first begin and last end of each super for length calculations.
    int maxSuperId = mrpSuperDataMgr->GetSize();

    mVecSuperBegin.resize( 0 );
    mVecSuperBegin.resize( maxSuperId + 1, 0 );
    
    mVecSuperEnd.resize( 0 );
    mVecSuperEnd.resize( maxSuperId + 1, 0 );

    mSumContigLengths.resize( 0 );
    mSumContigLengths.resize( maxSuperId + 1, 0 );

    for ( unsigned int clIdx = 0; clIdx < mVecContigLocations.size(); ++clIdx )
    {
        ContigLocation &contigLoc = mVecContigLocations[ clIdx ];
        int superId = contigLoc.GetSuper().GetId();
        
        mVecSuperBegin[ superId ] = min( mVecSuperBegin[ superId ], contigLoc.Begin() );
        mVecSuperEnd  [ superId ] = max( mVecSuperEnd[ superId ], contigLoc.End() );
        mSumContigLengths[ superId ] += contigLoc.Length();
    }

    mbDirty = false;

    // cout << "Reindexed ContigLocations in " << TimeSince( startTime ) << "." << endl;
}
