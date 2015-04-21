///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#include "assembly/ContigDataManager.h"

#include "assembly/ContigLocationManager.h"
#include "assembly/ContigSequenceManager.h"
#include "assembly/IdManager.h"
#include "assembly/ReadLocationManager.h"


ContigDataManager::ContigDataManager( ContigLocationManager *pCLM, 
                                      ReadLocationManager *pRLM )
    : mpIdMgr( new IdManager ),
      mpContigSequenceMgr( new ContigSequenceManager ),
      mpContigLocationMgr( pCLM ),
      mpReadLocationMgr( pRLM )
{
    ForceAssert( mpIdMgr );
    ForceAssert( mpContigSequenceMgr );
    ForceAssert( mpContigLocationMgr );
    ForceAssert( mpReadLocationMgr );
}

ContigDataManager::ContigDataManager( ContigLocationManager *pCLM, 
                                      ReadLocationManager *pRLM,
                                      ContigSequenceManager *pCSM, 
                                      IdManager *pIM )
    : mpIdMgr( pIM ),
      mpContigSequenceMgr( pCSM ),
      mpContigLocationMgr( pCLM ),
      mpReadLocationMgr( pRLM )
{
    ForceAssert( mpIdMgr );
    ForceAssert( mpContigSequenceMgr );
    ForceAssert( mpContigLocationMgr );
    ForceAssert( mpReadLocationMgr );
}

ContigDataManager::~ContigDataManager( )
{
    delete mpIdMgr;
    delete mpContigSequenceMgr;
}


void
ContigDataManager::SetLocationManager( ContigLocationManager *pCLM )
{
    mpContigLocationMgr = pCLM;
}


int
ContigDataManager::GetSize( ) const
{
    return mpIdMgr->GetMaxId() + 1;
}

Contig
ContigDataManager::NewContig( const basevector &bases, const qualvector &quals )
{
    int id = mpIdMgr->GetNewId();

    mpContigSequenceMgr->SetBases( id, bases );
    mpContigSequenceMgr->SetQuals( id, quals );

    mpContigSequenceMgr->Verify( id );

    return Contig( this, id );
}

bool
ContigDataManager::CheckId( const int id ) const
{
    ForceAssertGe( id, 0 );
    ForceAssertLe( id, mpIdMgr->GetMaxId() );
    return true;
}

Contig
ContigDataManager::GetContig( const int id )
{
  if ( id >= 0 && id <= mpIdMgr->GetMaxId() )
    return Contig( this, id );
  else
    return Contig();
}


int
ContigDataManager::GetLength( const int id ) const
{
    this->CheckId( id );

    return mpContigSequenceMgr->GetLength( id );
}

int
ContigDataManager::GetNumReadLocations( const int id ) const
{
    this->CheckId( id );

    return mpReadLocationMgr->GetNumReadLocationsInContig( id );
}


const basevector &
ContigDataManager::GetBases( const int id ) const
{
    this->CheckId( id );
    
    return mpContigSequenceMgr->GetBases( id );
}

const qualvector &
ContigDataManager::GetQuals( const int id ) const
{
    this->CheckId( id );
    
    return mpContigSequenceMgr->GetQuals( id );
}


void
ContigDataManager::SetBases( const int id, const basevector &newBases ) 
{
    this->CheckId( id );
    
    mpContigSequenceMgr->SetBases( id, newBases );
}

void
ContigDataManager::SetQuals( const int id, const qualvector &newQuals ) 
{
    this->CheckId( id );
    
    mpContigSequenceMgr->SetQuals( id, newQuals );
}


void
ContigDataManager::GetReadLocations( const int id, vec<ReadLocation> &vecRLs ) const
{
    mpReadLocationMgr->GetByContig( id, vecRLs );
}

void
ContigDataManager::GetContigLocations( const int id, vec<ContigLocation> &vecCLs ) const
{
    mpContigLocationMgr->GetByContig( id, vecCLs );
}

void
ContigDataManager::Reverse( const int id )
{
    this->CheckId( id );

    mpReadLocationMgr->Reverse( ContigToken( this, id ) );
    mpContigSequenceMgr->Reverse( id );
}

int
ContigDataManager::PurgeEmptyContigs()
{
  vec<ReadLocation> readLocs;
  vec<ContigLocation> contigLocs;

  int newContigId = 0, oldContigId;
  for ( oldContigId = 0; oldContigId < mpIdMgr->GetMaxId() + 1; ++oldContigId )
  {
    int contigLength = mpContigSequenceMgr->GetLength( oldContigId );
    int numReadsInContig = mpReadLocationMgr->GetNumReadLocationsInContig( oldContigId );

    if ( contigLength > 0 || numReadsInContig > 0 )
    {
      if ( contigLength == 0 )
        cout << "WARNING: C" << newContigId << " has reads but no bases." << endl;

      if ( numReadsInContig == 0 )
        cout << "WARNING: C" << newContigId << " has bases but no reads." << endl;

      if ( newContigId != oldContigId )
      {
        this->GetReadLocations( oldContigId, readLocs );
    
        for ( vec<ReadLocation>::iterator readLocIter = readLocs.begin();
              readLocIter != readLocs.end(); ++readLocIter )
        {
          mpReadLocationMgr->Remove( *readLocIter );
          mpReadLocationMgr->Add( ReadLocation( Contig( this, newContigId ),
                                                readLocIter->GetRead(),
                                                readLocIter->GetInterval(),
                                                readLocIter->GetOrientation() ) );
        }
        
        this->GetContigLocations( oldContigId, contigLocs );
        for ( vec<ContigLocation>::iterator contigLocIter = contigLocs.begin();
              contigLocIter != contigLocs.end(); ++contigLocIter )
        {
          mpContigLocationMgr->Remove( *contigLocIter );
          mpContigLocationMgr->Add( ContigLocation( contigLocIter->GetSuper(),
                                                    Contig( this, newContigId ),
                                                    contigLocIter->GetInterval(),
                                                    contigLocIter->GetOrientation() ) );
        }
        
        mpContigSequenceMgr->SetBases( newContigId, mpContigSequenceMgr->GetBases( oldContigId ) );
        mpContigSequenceMgr->SetQuals( newContigId, mpContigSequenceMgr->GetQuals( oldContigId ) );
      }

      ++newContigId;
    }
    else
    {
      this->GetContigLocations( oldContigId, contigLocs );
      for ( vec<ContigLocation>::iterator contigLocIter = contigLocs.begin();
            contigLocIter != contigLocs.end(); ++contigLocIter )
        mpContigLocationMgr->Remove( *contigLocIter );
    }      
  }

  delete mpIdMgr;
  mpIdMgr = new IdManager;
  mpIdMgr->UseId( newContigId-1 );

  return oldContigId - newContigId;
}

void
ContigDataManager::Write( const bool bOverwrite,
                          const String &strContigFastbFile,
                          const String &strContigQualbFile )
{
    mpContigSequenceMgr->Write( bOverwrite,
                                mpIdMgr->GetMaxId(),
                                strContigFastbFile,
                                strContigQualbFile );
}

void
ContigDataManager::SetPrefetchStrategy( PrefetchStrategy *pStrategy ) 
{
    mpContigSequenceMgr->SetPrefetchStrategy( pStrategy );
}
