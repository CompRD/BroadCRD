///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSEMBLY_READPAIRMANAGER
#define ASSEMBLY_READPAIRMANAGER

#include "String.h"
#include "Vec.h"

#include "assembly/ReadPairToken.h"

#include "ReadPairing.h"

class ReadDataManager;

class ReadPairManager
{
 public:
  ReadPairManager( ReadDataManager *&rpReadDataMgr, 
                   const String &strReadPairFile = "" );

  ReadPairToken Add( const ReadToken &read1, const ReadToken &read2, 
                     const int sep, const int stdev );

  int GetSize();

  ReadPairToken Get( const int pairId );
  ReadPairToken GetByRead( const int readId );
  
  ReadToken GetRead1( const int pairId ) const;
  ReadToken GetRead2( const int pairId ) const;
  
  int GetExpectedSep( const int pairId ) const;
  int GetExpectedStdDev( const int pairId ) const;
  
  void SetExpectedSep( const int pairId, const int sep );
  void SetExpectedStdDev( const int pairId, const int stdev );
  
  void IgnorePair( const int pairId );
  void UsePair( const int pairId );
  
  void Write( const bool bOverwrite,
              const String &strOutDir,
              const String &strReadPairsFile );
  
 private:
  void VerifyId( const int pairId ) const;

  int GetPairIdForRead( const int readId );
  
  bool mbLoaded;
  bool mbModified;
  
  String mStrReadPairFile;
  
  ReadDataManager * &mrpReadDataMgr;
  
  void Load();
  
  vec<read_pairing> mVecReadPairings;
  vec<int> mVecReadIdToPairMap;
};


// INLINE METHODS

inline
void 
ReadPairManager::VerifyId( const int pairId ) const
{
  AssertGe( pairId, 0 );
  AssertLt( pairId, (int) mVecReadPairings.size() );
  Assert( mVecReadPairings[ pairId ].Alive() );
}


inline
int
ReadPairManager::GetPairIdForRead( const int readId )
{
  if ( readId >= (int) mVecReadIdToPairMap.size() )
    return -1;

  return mVecReadIdToPairMap[ readId ];
}
  

inline
ReadPairToken
ReadPairManager::Get( const int pairId ) 
{
  if ( ! mbLoaded )
    this->Load();

  if ( pairId < 0 || mVecReadPairings[ pairId ].Dead() )
    return ReadPairToken();
  else
    return ReadPairToken( this, pairId );
}


inline
ReadPairToken
ReadPairManager::GetByRead( const int readId ) 
{
  if ( ! mbLoaded )
    this->Load();

  int pairId = GetPairIdForRead( readId );
  if ( pairId < 0 || mVecReadPairings[ pairId ].Dead() )
    return ReadPairToken();
  else
    return ReadPairToken( this, pairId );
}


inline
ReadToken
ReadPairManager::GetRead1( const int pairId ) const
{
  VerifyId( pairId );
  return ReadToken( mrpReadDataMgr, mVecReadPairings[ pairId ].id1 );
}

inline
ReadToken
ReadPairManager::GetRead2( const int pairId ) const
{
  VerifyId( pairId );
  return ReadToken( mrpReadDataMgr, mVecReadPairings[ pairId ].id2 );
}

inline
int 
ReadPairManager::GetExpectedSep( const int pairId ) const
{
  VerifyId( pairId );
  return mVecReadPairings[ pairId ].sep;
}

inline
int
ReadPairManager::GetExpectedStdDev( const int pairId ) const
{
  VerifyId( pairId );
  return mVecReadPairings[ pairId ].sd;
}


inline
void
ReadPairManager::SetExpectedSep( const int pairId, const int sep ) 
{
  VerifyId( pairId );
  mVecReadPairings[ pairId ].sep = sep;
}

inline
void
ReadPairManager::SetExpectedStdDev( const int pairId, const int stdev )
{
  VerifyId( pairId );
  mVecReadPairings[ pairId ].sd = stdev;
}


#endif
