///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "assembly/SuperDataManager.h"

#include "assembly/IdManager.h"
#include "assembly/ContigLocationManager.h"

SuperDataManager::SuperDataManager( ContigLocationManager *pCLM )
    : mpIdMgr( new IdManager ),
      mpContigLocationMgr( pCLM )
{
    ForceAssert( mpIdMgr );
    ForceAssert( mpContigLocationMgr );
}

SuperDataManager::SuperDataManager( ContigLocationManager *pCLM,
                                    IdManager *pIM )
    : mpIdMgr( pIM ),
      mpContigLocationMgr( pCLM )
{
    ForceAssert( mpIdMgr );
    ForceAssert( mpContigLocationMgr );
}

SuperDataManager::~SuperDataManager( )
{
    delete mpIdMgr;
}


int
SuperDataManager::GetSize( ) const
{
    return mpIdMgr->GetMaxId() + 1;
}

Super
SuperDataManager::NewSuper( ) 
{
    int id = mpIdMgr->GetNewId();

    return Super( this, id );
}

Super
SuperDataManager::GetSuper( const int id )
{
  if ( id >= 0 && id <= mpIdMgr->GetMaxId() )
    return Super( this, id );
  else
    return Super();
}


int
SuperDataManager::GetNumContigLocations( const int id ) const
{
    return mpContigLocationMgr->GetNumContigLocationsInSuper( id );
}

void
SuperDataManager::GetContigLocations( const int id, vec<ContigLocation> &vecCLs ) const
{
    mpContigLocationMgr->GetBySuper( id, vecCLs );
}


int 
SuperDataManager::GetLength( const int id ) const
{
    return mpContigLocationMgr->GetSuperLength( id );
}

longlong
SuperDataManager::GetSumOfContigLengths( const int id ) const
{
    return mpContigLocationMgr->GetSumOfContigLengths( id );
}

void
SuperDataManager::Reverse( const int id )
{
    ForceAssertGe( id, 0 );
    ForceAssertLe( id, mpIdMgr->GetMaxId() );

    mpContigLocationMgr->Reverse( SuperToken( this, id ) );
}
