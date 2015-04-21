///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "assembly/ContigSequenceManager.h"

#include "assembly/PrefetchStrategy.h"

ContigSequenceManager::ContigSequenceManager( )
    : mpBasesMgr( new BasesManager ),
      mpQualsMgr( new QualsManager ),
      mpPrefetchStrategy( new NullPrefetchStrategy )
{
}

ContigSequenceManager::ContigSequenceManager( const String &strContigFastbFile,
                                              const String &strContigQualbFile )
    : mpBasesMgr( new BasesManager( strContigFastbFile ) ),
      mpQualsMgr( new QualsManager( strContigQualbFile ) ),
      mpPrefetchStrategy( new NullPrefetchStrategy )
{
}

ContigSequenceManager::~ContigSequenceManager( )
{
    delete mpBasesMgr;
    delete mpQualsMgr;
    delete mpPrefetchStrategy;
}

void
ContigSequenceManager::LoadBases( const int id ) const
{
    mpBasesMgr->LoadData( id );
}

void
ContigSequenceManager::LoadQuals( const int id ) const
{
    mpQualsMgr->LoadData( id );
}


void
ContigSequenceManager::LoadBases( const vec<int> vecIds ) const
{
    mpBasesMgr->LoadData( vecIds );
}

void
ContigSequenceManager::LoadQuals( const vec<int> vecIds ) const
{
    mpQualsMgr->LoadData( vecIds );
}


void
ContigSequenceManager::LoadAllBases( ) const
{
    mpBasesMgr->LoadAllData();
}

void
ContigSequenceManager::LoadAllQuals( ) const
{
    mpQualsMgr->LoadAllData( );
}


const basevector &
ContigSequenceManager::GetBases( const int id ) const
{
    if ( ! mpBasesMgr->IsLoaded( id ) )
    {
        if ( mpPrefetchStrategy->GetAll() )
            mpBasesMgr->LoadAllData();
        else
        {
            vec<int> vecIds;
            mpPrefetchStrategy->GetIdsToLoad( id, vecIds );
            mpBasesMgr->LoadData( vecIds );
        }
    }

    return mpBasesMgr->GetData( id );
}

const qualvector &
ContigSequenceManager::GetQuals( const int id ) const
{
    if ( ! mpQualsMgr->IsLoaded( id ) )
    {
        if ( mpPrefetchStrategy->GetAll() )
            mpQualsMgr->LoadAllData();
        else
        {
            vec<int> vecIds;
            mpPrefetchStrategy->GetIdsToLoad( id, vecIds );
            mpQualsMgr->LoadData( vecIds );
        }
    }

    return mpQualsMgr->GetData( id );
}


void
ContigSequenceManager::SetBases( const int id, const basevector &bases )
{
    mpBasesMgr->SetData( id, bases );
}

void
ContigSequenceManager::SetQuals( const int id, const qualvector &quals )
{
    mpQualsMgr->SetData( id, quals );
}


bool
ContigSequenceManager::Verify( const int id ) const
{
    ForceAssertEq( GetBases(id).size(), GetQuals(id).size() );

    return true;
}

int
ContigSequenceManager::GetLength( const int id ) const
{
    return GetBases( id ).size();
}

void
ContigSequenceManager::Reverse( const int id )
{
    basevector bases( GetBases( id ) );
    bases.ReverseComplement();

    SetBases( id, bases );

    qualvector quals( GetQuals( id ) );
    quals.ReverseMe();

    SetQuals( id, quals );
}

void
ContigSequenceManager::Write( const bool bOverwrite,
                              const int lastId,
                              const String &strContigFastbFile,
                              const String &strContigQualbFile )
{
    mpBasesMgr->Write( bOverwrite, strContigFastbFile, lastId );
    mpQualsMgr->Write( bOverwrite, strContigQualbFile, lastId );
}

void
ContigSequenceManager::SetPrefetchStrategy( PrefetchStrategy *pPrefetchStrategy )
{
    delete mpPrefetchStrategy;
    mpPrefetchStrategy = pPrefetchStrategy->Clone();
}
