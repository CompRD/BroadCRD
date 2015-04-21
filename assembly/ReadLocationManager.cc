///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "ReadLocation.h"
#include "STLExtensions.h"

#include "assembly/ReadLocationManager.h"
#include "assembly/ContigDataManager.h"
#include "assembly/ReadDataManager.h"

const int k_readIsUnplaced = -1;
const int k_readIsMultiPlaced = -2;

ReadLocationManager::ReadLocationManager( ContigDataManager * &rpContigDataMgr,
                                          ReadDataManager * &rpReadDataMgr,
                                          const String &strReadLocationFile )
    : mbLoaded( strReadLocationFile.empty() ),
      mbModified( false ),
      mStrReadLocationFile( strReadLocationFile ),
      mrpContigDataMgr( rpContigDataMgr ),
      mrpReadDataMgr( rpReadDataMgr ),
      mbByContigIsOutOfDate( true ),
      mbByReadIsOutOfDate( true )
{
}


void
ReadLocationManager::Load()
{
    if ( mbLoaded )
        return;

    // cout << "Filling out ReadLocationManager." << endl;

    READ( mStrReadLocationFile, vec<read_location>, vecReadLocations );

    mVecReadLocations.resize( 0 );
    mVecReadLocations.reserve( vecReadLocations.size() );
    
    for ( int locIdx = 0; locIdx < (int) vecReadLocations.size(); ++locIdx )
    {
        read_location &theLoc = vecReadLocations[ locIdx ];
        
        ContigToken theContig = ContigToken( mrpContigDataMgr, theLoc.Contig() );
        ReadToken   theRead   = ReadToken( mrpReadDataMgr, theLoc.ReadId() );
        int         begin     = theLoc.StartOnContig();
        int         end       = theLoc.StopOnContig() + 1;
        Orientation theOrient = 
            ( theLoc.OrientationOnContig() == ForwardOr ? orient_FW : orient_RC );
        
        mVecReadLocations.push_back( ReadLocation( theContig, theRead,
                                                   Interval( begin, end ), theOrient ) );
    }
    
    mbLoaded = true;
    mbByContigIsOutOfDate = true;
    mbByReadIsOutOfDate = true;
}



void
ReadLocationManager::Add( const ReadLocation &theReadLoc )
{
    if ( ! mbLoaded ) Load();

    mbModified = true;
    mVecReadLocations.push_back( theReadLoc );

    int newLocIdx = mVecReadLocations.size() - 1;

    int readId = theReadLoc.GetRead().GetId();
    if ( readId < 0 ) 
      return;

    while ( (int) mIndicesByRead.size() < readId + 1 )
      mIndicesByRead.push_back( k_readIsUnplaced );

    if ( mIndicesByRead[ readId ] == k_readIsUnplaced )
      mIndicesByRead[ readId ] = newLocIdx;

    else 
    {
      if ( mIndicesByRead[ readId ] >= 0 )
      {
        mMapMultiPlacedReads.insert( make_pair( readId, mIndicesByRead[ readId ] ) );
        mIndicesByRead[ readId ] = k_readIsMultiPlaced;
      }

      ForceAssertEq( mIndicesByRead[ readId ], k_readIsMultiPlaced );
      mMapMultiPlacedReads.insert( make_pair( readId, newLocIdx ) );
    }

    int contigId = theReadLoc.GetContig().GetId();
    if ( contigId < 0 ) 
      return;
    
    while ( (int) mIndicesByContig.size() < contigId + 1 )
      mIndicesByContig.push_back( vec<int>() );

    mIndicesByContig[ contigId ].push_back( newLocIdx );
}


bool
ReadLocationManager::Remove( const ReadLocation &theReadLoc )
{
    if ( ! mbLoaded ) Load();
    if ( mbByReadIsOutOfDate ) UpdateByRead();

    int locIdx = -1;

    int theReadId = theReadLoc.GetRead().GetId();

    // Find the index of the ReadLocation by checking the read-to-loc map.

    locIdx = mIndicesByRead[ theReadId ];

    if ( locIdx >= 0 )
      mIndicesByRead[ theReadId ] = k_readIsUnplaced;

    else if ( locIdx == k_readIsUnplaced )
      return false;

    else
    {
      ForceAssertEq( locIdx, k_readIsMultiPlaced );
      
      pair<multimap<int,int>::iterator,multimap<int,int>::iterator> mapRange;
      mapRange = mMapMultiPlacedReads.equal_range( theReadId );

      for ( ; mapRange.first != mapRange.second; ++mapRange.first )
      {
        locIdx = mapRange.first->second;
        if ( mVecReadLocations[ locIdx ] == theReadLoc )
          break;
      }

      if ( mapRange.first == mapRange.second )
        return false;

      else
        mMapMultiPlacedReads.erase( mapRange.first );

      // Is the read singly-placed now?
      mapRange = mMapMultiPlacedReads.equal_range( theReadId );

      if ( distance( mapRange.first, mapRange.second ) == 1 )
      {
        mIndicesByRead[ theReadId ] = mapRange.first->second;
        mMapMultiPlacedReads.erase( mapRange.first, mapRange.second );
      }
    }

    if ( mbByContigIsOutOfDate ) UpdateByContig();

    Contig theContig = theReadLoc.GetContig();

    vec<int> &vecIndicesFromContig = mIndicesByContig[ theContig.GetId() ];

    vec<int>::iterator idxIter = vecIndicesFromContig.begin(); 
    while ( idxIter != vecIndicesFromContig.end() )
    {
      if ( *idxIter == locIdx )
      {
        vecIndicesFromContig.erase( idxIter );
        mbModified = true;
      }
      else
        ++idxIter;
    }

    mVecReadLocations[ locIdx ] = ReadLocation();


    return true;
}


void
ReadLocationManager::ClearContig( const int contigId )
{
    if ( mbByContigIsOutOfDate ) UpdateByContig();

    // If the contig is empty and has been added since we last called
    // UpdateByContig(), the index vector will be too small.  Check
    // that the id is legal and return.
    if ( contigId >= (int) mIndicesByContig.size() )
    {
      ForceAssertLt( contigId, mrpContigDataMgr->GetSize() );
      return;
    }

    vec<int> &vecIndices = mIndicesByContig[ contigId ];

    for ( vec<int>::reverse_iterator idxIter = vecIndices.rbegin(); 
          idxIter != vecIndices.rend(); ++idxIter )
        this->Remove( mVecReadLocations[ *idxIter ] );
}


void
ReadLocationManager::Reverse( const ContigToken theContig )
{
    if ( ! mbLoaded ) Load();
    if ( mbByContigIsOutOfDate ) UpdateByContig();

    // If the contig is empty and has been added since we last called
    // UpdateByContig(), the index vector will be too small.  Check
    // that the id is legal and return.

    const int contigId = theContig.GetId();

    if ( contigId >= (int) mIndicesByContig.size() )
    {
      ForceAssertLt( contigId, mrpContigDataMgr->GetSize() );
      return;
    }

    mbModified = true;

    vec<int> &indicesOfContig = mIndicesByContig[ contigId ];

    int length = Contig(theContig).GetLength();

    for ( vec<int>::iterator idxIter = indicesOfContig.begin();
          idxIter != indicesOfContig.end(); ++idxIter )
    {
        ReadLocation &theLoc = mVecReadLocations[ *idxIter ];

        theLoc = ReadLocation( theLoc.GetContig(), theLoc.GetRead(),
                               Interval( length - theLoc.End(), length - theLoc.Begin() ),
                               Flip( theLoc.GetOrientation() ) );
    }
}

    
void
ReadLocationManager::GetByContig( const int contigId,
                                  vec<ReadLocation> &vecReadLocs )
{
    if ( ! mbLoaded ) Load();
    if ( mbByContigIsOutOfDate ) UpdateByContig();

    vecReadLocs.clear();

    // If the contig is empty and has been added since we last called
    // UpdateByContig(), the index vector will be too small.  Check
    // that the id is legal and return.
    if ( contigId >= (int) mIndicesByContig.size() )
    {
      ForceAssertLt( contigId, mrpContigDataMgr->GetSize() );
      return;
    }

    vec<int> &indicesOfContig = mIndicesByContig[ contigId ];

    vecReadLocs.reserve( indicesOfContig.size() );

    for ( vec<int>::iterator idxIter = indicesOfContig.begin();
          idxIter != indicesOfContig.end(); ++idxIter )
        vecReadLocs.push_back( mVecReadLocations[ *idxIter ] );
}

void
ReadLocationManager::GetByRead( const int readId, 
                                vec<ReadLocation> &vecReadLocs )
{
    if ( ! mbLoaded ) Load();
    if ( mbByReadIsOutOfDate ) UpdateByRead();

    vecReadLocs.clear();

    // If the read is unplaced and has been added since we last called
    // UpdateByRead(), the index vector will be too small.  Check
    // that the id is legal and return.
    if ( readId >= (int) mIndicesByRead.size() )
    {
      ForceAssertLt( readId, mrpReadDataMgr->GetSize() );
      return;
    }

    int rlIdx = mIndicesByRead[ readId ];

    // rlIdx is either -1 (unplaced), -2 (multiply placed), or the
    // locIdx of the read's one placement.
    
    if ( rlIdx >= 0 )
    {
      vecReadLocs.push_back( mVecReadLocations[ rlIdx ] );
      return;
    }

    else if ( rlIdx == k_readIsUnplaced )
    {
      return;
    }

    else // rlIdx == k_readIsMultiPlaced
    {
      ForceAssertEq( rlIdx, k_readIsMultiPlaced );
      
      pair<multimap<int,int>::iterator,multimap<int,int>::iterator> range;
      
      range = mMapMultiPlacedReads.equal_range( readId );

      for ( ; range.first != range.second; ++range.first )
        vecReadLocs.push_back( mVecReadLocations[ range.first->second ] );
      
      return;
    }
}


int
ReadLocationManager::GetNumReadLocationsOfRead( const int readId )
{
    if ( ! mbLoaded ) Load();
    if ( mbByReadIsOutOfDate ) UpdateByRead();

    int rlIdx = mIndicesByRead[ readId ];

    if ( rlIdx == k_readIsUnplaced )
      return 0;

    else if ( rlIdx >= 0 )
      return 1;

    else // rlIdx == k_readIsMultiPlaced
    {
      ForceAssertEq( rlIdx, k_readIsMultiPlaced );

      pair<multimap<int,int>::iterator,multimap<int,int>::iterator> range;
      range = mMapMultiPlacedReads.equal_range( readId );

      return distance( range.first, range.second );
    }
}


int
ReadLocationManager::GetNumReadLocationsInContig( const int contigId )
{
    if ( ! mbLoaded ) Load();
    if ( mbByContigIsOutOfDate ) UpdateByContig();

    // If the contig is empty and has been added since we last called
    // UpdateByContig(), the index vector will be too small.  Check
    // that the id is legal and return.
    if ( contigId >= (int) mIndicesByContig.size() )
    {
      ForceAssertLt( contigId, mrpContigDataMgr->GetSize() );
      return 0;
    }

    const vec<int> &indicesOfContig = mIndicesByContig[ contigId ];
    return indicesOfContig.size();
}


void 
ReadLocationManager::UpdateByContig()
{
    if ( ! mbLoaded ) Load();

    // cout << "Indexing ReadLocationManager by contig." << endl;

    for ( unsigned int idx = 0; idx < mIndicesByContig.size(); ++idx )
        mIndicesByContig[ idx ].resize( 0 );

    mIndicesByContig.resize( mrpContigDataMgr->GetSize() );

    for ( unsigned int rlIdx = 0; rlIdx < mVecReadLocations.size(); ++rlIdx )
    {
        ReadLocation &theReadLoc = mVecReadLocations[ rlIdx ];

        if ( ! theReadLoc.GetRead().IsValid() )
          continue;

        int contigId = theReadLoc.GetContig().GetId();

        mIndicesByContig[ contigId ].push_back( rlIdx );
    }

    mbByContigIsOutOfDate = false;
}

void 
ReadLocationManager::UpdateByRead()
{
    if ( ! mbLoaded ) Load();

    // cout << "Indexing ReadLocationManager by read." << endl;

    // the entry for a given read id in mIndicesByRead is either -1
    // (unplaced), -2 (multiply placed), or the locIdx of the read's
    // one placement.
    
    mIndicesByRead.resize( mrpReadDataMgr->GetSize() );
    fill( mIndicesByRead.begin(), mIndicesByRead.end(), k_readIsUnplaced );

    mMapMultiPlacedReads.clear();

    for ( unsigned int rlIdx = 0; rlIdx < mVecReadLocations.size(); ++rlIdx )
    {
        ReadLocation &theReadLoc = mVecReadLocations[ rlIdx ];

        if ( ! theReadLoc.GetRead().IsValid() )
          continue;

        int readId = theReadLoc.GetRead().GetId();

        if ( mIndicesByRead[ readId ] == k_readIsUnplaced )
          mIndicesByRead[ readId ] = rlIdx;

        else if ( mIndicesByRead[ readId ] >= 0 )
        {
          mMapMultiPlacedReads.insert( make_pair( readId, mIndicesByRead[ readId ] ) );
          mMapMultiPlacedReads.insert( make_pair( readId, rlIdx ) );
          mIndicesByRead[ readId ] = k_readIsMultiPlaced;
        }

        else // mIndicesByRead[ readId ] == k_readIsMultiPlaced
        {
          ForceAssertEq( mIndicesByRead[ readId ], k_readIsMultiPlaced );
          mMapMultiPlacedReads.insert( make_pair( readId, rlIdx ) );
        }
    }

    mbByReadIsOutOfDate = false;
}

    
void
ReadLocationManager::Write( const bool bOverwrite,
                            const String &strReadLocationFile )
{
    if ( ! mbModified && 
         IsRegularFile(strReadLocationFile) && 
         IsRegularFile(mStrReadLocationFile) &&
         RealPath(strReadLocationFile) == RealPath(mStrReadLocationFile) )
        return;

    ForceAssert( bOverwrite || ! IsRegularFile( strReadLocationFile ) );
    
    this->Load();

    vec<read_location> vecOrigReadLocs;
    vecOrigReadLocs.reserve( mVecReadLocations.size() );

    for ( int locIdx = 0; locIdx < (int) mVecReadLocations.size(); ++locIdx )
    {
        ReadLocation &theLoc = mVecReadLocations[ locIdx ];

        if ( ! theLoc.GetRead().IsValid() )
          continue;

        read_location origLoc( theLoc.GetRead().GetId(),
                               theLoc.Length(),
                               theLoc.GetContig().GetId(),
                               theLoc.Begin(),
                               ( theLoc.GetOrientation() == orient_FW ? ForwardOr : ReverseOr ),
                               Contig(theLoc.GetContig()).GetLength() );
        
        vecOrigReadLocs.push_back( origLoc );
    }

    sort( vecOrigReadLocs.begin(), vecOrigReadLocs.end() );

    int numContigs = mrpContigDataMgr->GetSize();
    int numReads   = mrpReadDataMgr->GetSize();

    if ( vecOrigReadLocs.empty() )
    {
      cout << "WARNING: there are no reads placed in this assembly." << endl;

      // Set numContigs and numReads to special values that prevent WriteLocs() from asserting.
      numContigs = -1;
      numReads = -1;
    }

    WriteLocs( strReadLocationFile, vecOrigReadLocs, 
               numContigs, numReads );
}
