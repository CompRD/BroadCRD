///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "assembly/NameManager.h"
#include "assembly/PrefetchStrategy.h"

NameManager::NameManager()
    : mpStringsMgr( new StringsManager ),
      mbDirty( false ),
      mpPrefetchStrategy( new NullPrefetchStrategy )
{
}

NameManager::NameManager( const String &strNamesFile )
    : mpStringsMgr( new StringsManager( strNamesFile ) ),
      mbDirty( true ),
      mpPrefetchStrategy( new NullPrefetchStrategy )
{
}

NameManager::~NameManager( )
{
    delete mpStringsMgr;
    delete mpPrefetchStrategy;
}

String
NameManager::GetName( const int id ) const
{
    if ( ! mpStringsMgr->IsLoaded( id ) )
    {
        if ( mpPrefetchStrategy->GetAll() )
            mpStringsMgr->LoadAllData();
        else
        {
            vec<int> vecIds;
            mpPrefetchStrategy->GetIdsToLoad( id, vecIds );
            mpStringsMgr->LoadData( vecIds );
        }
    }

    return mpStringsMgr->GetData( id );
}


void
NameManager::SetName( const int id, const String &name )
{
    mpStringsMgr->SetData( id, name );
    mMapNameToId.insert( make_pair( mpStringsMgr->GetData( id ), id ) );
}

int
NameManager::FindId( const String &name ) const
{
    if ( IsDirty() ) Clean();

    map<String,int>::const_iterator map_iter = mMapNameToId.find( name );

    if ( map_iter == mMapNameToId.end() )
        return -1;

    return map_iter->second;
}
                               
bool
NameManager::Verify( const int id, AddBehavior b ) const
{
    const String &name = mpStringsMgr->GetData( id );
    ForceAssert( ! name.empty() );

    int foundId = FindId( name );

    if ( b == ASSERT_IF_DIFFERENT )
      ForceAssertEq( id, foundId );

    return true;
}

void
NameManager::Write( const bool bOverwrite,
                        const String &strNamesFile )
{
    mpStringsMgr->Write( bOverwrite,
                         strNamesFile );
}

bool 
NameManager::IsDirty( ) const
{
    return ( mbDirty );
}

void
NameManager::Clean( ) const
{
    mpStringsMgr->LoadAllData();
    int numNames = mpStringsMgr->Size();
    for ( int id = 0; id < numNames; ++id )
        mMapNameToId.insert( make_pair( mpStringsMgr->GetData( id ), id ) );
    
    mbDirty = false;
}

void
NameManager::SetPrefetchStrategy( PrefetchStrategy *pPrefetchStrategy )
{
    delete mpPrefetchStrategy;
    mpPrefetchStrategy = pPrefetchStrategy->Clone();
}
