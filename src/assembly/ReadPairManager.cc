///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#include "system/Assert.h"
#include "ReadPairing.h"

#include "assembly/ReadDataManager.h"
#include "assembly/ReadPairManager.h"

ReadPairManager::ReadPairManager( ReadDataManager * &rpReadDataMgr,
                                  const String &strReadPairFile )
    : mbLoaded( strReadPairFile.empty() ),
      mbModified( false ),
      mStrReadPairFile( strReadPairFile ),
      mrpReadDataMgr( rpReadDataMgr )
{
}   

int
ReadPairManager::GetSize()
{
  if ( ! mbLoaded )
    this->Load();

  return mVecReadPairings.size();
}

ReadPairToken
ReadPairManager::Add( const ReadToken &read1, const ReadToken &read2,
                      const int sep, const int stdev )
{
  ForceAssert( read1.GetMgrPtr() == mrpReadDataMgr );
  ForceAssert( read2.GetMgrPtr() == mrpReadDataMgr );
  ForceAssertEq( GetPairIdForRead( read1.GetId() ), -1 );
  ForceAssertEq( GetPairIdForRead( read2.GetId() ), -1 );

  read_pairing newPair;
  newPair.id1 = read1.GetId();
  newPair.id2 = read2.GetId();
  newPair.sep = sep;
  newPair.sd  = stdev;

  int newPairId = mVecReadPairings.size();
  mVecReadPairings.push_back( newPair );

  int minMapSize = max( read1.GetId(), read2.GetId() ) + 1;
  if ( (int) mVecReadIdToPairMap.size() < minMapSize )
    mVecReadIdToPairMap.resize( minMapSize, -1 );

  mVecReadIdToPairMap[ read1.GetId() ] = newPairId;
  mVecReadIdToPairMap[ read2.GetId() ] = newPairId;

  return ReadPairToken( this, newPairId );
}

void
ReadPairManager::Load()
{
    //    cout << "Filling out ReadPairManager... " << flush;

    mVecReadIdToPairMap.resize( 0 );
    mVecReadIdToPairMap.resize( mrpReadDataMgr->GetSize(), -1 );
    
    ReadPairsFile( mStrReadPairFile, mVecReadPairings );

    for ( int pairId = 0; pairId < (int) mVecReadPairings.size(); ++pairId )
    {
        read_pairing &thePair = mVecReadPairings[ pairId ];
        
        if ( thePair.Dead() )
          continue;

        mVecReadIdToPairMap[ thePair.id1 ] = pairId;
        mVecReadIdToPairMap[ thePair.id2 ] = pairId;
    }

    mbLoaded = true;

    //    cout << "done." << endl;
}    


void
ReadPairManager::IgnorePair( const int pairId ) 
{
  VerifyId( pairId );

  read_pairing &thePair = mVecReadPairings[ pairId ];
  
  if ( thePair.Dead() )
    return;

  mVecReadIdToPairMap[ thePair.id1 ] = -1;
  mVecReadIdToPairMap[ thePair.id2 ] = -1;
}


void
ReadPairManager::UsePair( const int pairId ) 
{
  VerifyId( pairId );

  read_pairing &thePair = mVecReadPairings[ pairId ];
  
  if ( thePair.Dead() )
    return;

  mVecReadIdToPairMap[ thePair.id1 ] = pairId;
  mVecReadIdToPairMap[ thePair.id2 ] = pairId;
}


void
ReadPairManager::Write( const bool bOverwrite,
                        const String &strOutDir,
                        const String &strOutFile )
{
  ForceAssert( strOutFile.Contains( strOutDir, 0 ) );

  if ( ! mbModified &&          
       IsRegularFile(strOutFile) && 
       IsRegularFile(mStrReadPairFile) &&
       RealPath(strOutFile) == RealPath(mStrReadPairFile) )
    return;
  
  ForceAssert( bOverwrite || ! IsRegularFile( strOutFile ) );
  
  if ( ! mbLoaded )
    this->Load();
  
  WritePairs( strOutDir, mVecReadPairings, mrpReadDataMgr->GetSize() );
}
