///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "assembly/ReadDataManager.h"

#include "assembly/BinaryVecDataManager.h"
#include "assembly/IdManager.h"
#include "assembly/NameManager.h"
#include "assembly/ReadSequenceManager.h"
#include "assembly/ReadLocationManager.h"
#include "assembly/ReadPairManager.h"


ReadDataManager::ReadDataManager( ReadLocationManager *pRLM, 
                                  ReadPairManager *pRPM )
    : mpIdMgr( new IdManager ),
      mpReadNameMgr( new NameManager ),
      mpReadSequenceMgr( new ReadSequenceManager ),
      mpRepetitiveFlagMgr( new BinaryVecDataManager<Bool> ),
      mpReadLocationMgr( pRLM ),
      mpReadPairMgr( pRPM )
{
    ForceAssert( mpIdMgr );
    ForceAssert( mpReadNameMgr );
    ForceAssert( mpReadSequenceMgr );
    ForceAssert( mpRepetitiveFlagMgr );
    ForceAssert( mpReadLocationMgr );
    ForceAssert( mpReadPairMgr );
}

ReadDataManager::ReadDataManager( ReadLocationManager *pRLM, 
                                  ReadPairManager *pRPM,
                                  NameManager *pRNM, 
                                  ReadSequenceManager *pRSM,
                                  BinaryVecDataManager<Bool> *pRFM,
                                  IdManager *pIM )
    : mpIdMgr( pIM ),
      mpReadNameMgr( pRNM ),
      mpReadSequenceMgr( pRSM ),
      mpRepetitiveFlagMgr( pRFM ),
      mpReadLocationMgr( pRLM ),
      mpReadPairMgr( pRPM )
{
    ForceAssert( mpIdMgr );
    ForceAssert( mpReadNameMgr );
    ForceAssert( mpReadSequenceMgr );
    ForceAssert( mpRepetitiveFlagMgr );
    ForceAssert( mpReadLocationMgr );
    ForceAssert( mpReadPairMgr );
}

ReadDataManager::~ReadDataManager( )
{
    delete mpIdMgr;
    delete mpReadNameMgr;
    delete mpReadSequenceMgr;
    delete mpRepetitiveFlagMgr;
}


void
ReadDataManager::SetLocationManager( ReadLocationManager *pRLM )
{
    mpReadLocationMgr = pRLM;
}


int
ReadDataManager::GetSize( ) const
{
    return mpIdMgr->GetMaxId() + 1;
}


Read
ReadDataManager::NewRead( const String &name,
                          const basevector &bases,
                          const qualvector &quals,
                          const pair<int,int> trims,
                          const CompressedSequence &untrimmedBases,
			  AddBehavior behavior )
{
    ForceAssert( ! name.empty() );

    Read existingRead = this->GetReadNamed( name );

    if ( existingRead.IsValid() && behavior != REPLACE_EXISTING )
    {
      if (   ( behavior == ASSERT_IF_DIFFERENT ) && 
	     ( ! ( existingRead.GetBases() == bases ) ) ||
             ( ! ( existingRead.GetQuals() == quals ) ) ||
             ( ! ( existingRead.GetTrims() == trims ) ) ||
             ( ! ( existingRead.GetUntrimmedBases() == untrimmedBases ) ) )
        {
            cerr << "Attempted to add read with same name but different data." << endl;
            PRINT2_TO( cerr, existingRead.GetBases().size(), bases.size() );
            // if ( existingRead.GetBases().size() == bases.size() )
            //     for ( unsigned int i = 0; i < bases.size(); ++i )
            //         if ( existingRead.GetBases()(i) != bases(i) )
            //             PRINT3_TO( cerr, 
            //                        as_base( existingRead.GetBases()(i) ),
            //                        as_base( bases(i) ), 
            //                        i );

            PRINT2_TO( cerr, existingRead.GetQuals().size(), quals.size() );
            // if ( existingRead.GetQuals().size() == quals.size() )
            //     for ( unsigned int i = 0; i < quals.size(); ++i )
            //         if ( existingRead.GetQuals()[i] != quals[i] )
            //             PRINT3_TO( cerr,
            //                        (int) existingRead.GetQuals()[i],
            //                        (int) quals[i],
            //                        i );

            PRINT2_TO( cerr, existingRead.GetTrims().first, trims.first );
            PRINT2_TO( cerr, existingRead.GetTrims().second, trims.second );

            PRINT2_TO( cerr, existingRead.GetUntrimmedBases().size(), untrimmedBases.size() );

	    
            ForceAssert( 0 == 1 );
            return Read();
        }
        else
            return existingRead;
    }
    else
    {
        int id = mpIdMgr->GetNewId( );

        mpReadNameMgr->SetName( id, name );
        mpReadSequenceMgr->SetBases( id, bases );
        mpReadSequenceMgr->SetQuals( id, quals );
        mpReadSequenceMgr->SetTrims( id, trims );
        mpReadSequenceMgr->SetUntrimmedBases( id, untrimmedBases );
        mpRepetitiveFlagMgr->SetData( id, False );

        mpReadNameMgr->Verify( id, behavior );
        mpReadSequenceMgr->Verify( id );

        return Read( this, id );
    }
}


bool
ReadDataManager::CheckId( const int id ) const
{
    ForceAssertGe( id, 0 );
    ForceAssertLe( id, mpIdMgr->GetMaxId() );
    return true;
}

Read
ReadDataManager::GetRead( const int id )
{
  if ( id >= 0 && id <= mpIdMgr->GetMaxId() )
    return Read( this, id );
  else
    return Read();
}


Read
ReadDataManager::GetReadNamed( const String &name ) 
{
    int id = mpReadNameMgr->FindId( name );

    if ( id == -1 )
        return Read();

    return this->GetRead( id );
}
    

String
ReadDataManager::GetName( const int id ) const
{
    CheckId( id );
    
    return mpReadNameMgr->GetName( id );
}

int
ReadDataManager::GetLength( const int id ) const
{
    CheckId( id );
    
    return mpReadSequenceMgr->GetLength( id );
}

const basevector &
ReadDataManager::GetBases( const int id ) const
{
    CheckId( id );

    return mpReadSequenceMgr->GetBases( id );
}
    
const qualvector &
ReadDataManager::GetQuals( const int id ) const
{
    CheckId( id );

    return mpReadSequenceMgr->GetQuals( id );
}
    
const CompressedSequence &
ReadDataManager::GetUntrimmedBases( const int id ) const
{
    CheckId( id );

    return mpReadSequenceMgr->GetUntrimmedBases( id );
}
    
pair<int,int>
ReadDataManager::GetTrims( const int id ) const
{
    CheckId( id );

    return mpReadSequenceMgr->GetTrims( id );
}
    
Bool
ReadDataManager::IsRepetitive( const int id ) const
{
    CheckId( id );
  
    if ( !mpRepetitiveFlagMgr->fileExists() )
      return False;

    return mpRepetitiveFlagMgr->GetData( id );
}

int
ReadDataManager::GetNumLocations( const int id ) const
{
    CheckId( id );

    return mpReadLocationMgr->GetNumReadLocationsOfRead( id );
}

void
ReadDataManager::GetLocations( const int id, vec<ReadLocation> &vecRLs ) const
{
    CheckId( id );

    mpReadLocationMgr->GetByRead( id, vecRLs );
}

ReadPair
ReadDataManager::GetPair( const int id ) const
{
    CheckId( id );
    
    return mpReadPairMgr->GetByRead( id );
}

void
ReadDataManager::Write( const bool bOverwrite,
                        const String &strReadNamesFile,
                        const String &strReadLengthsFile,
                        const String &strReadFastbFile,
                        const String &strReadQualbFile,
                        const String &strReadTrimsFile,
                        const String &strReadFastnFile )
{
    mpReadNameMgr->Write( bOverwrite, 
                          strReadNamesFile );

    mpReadSequenceMgr->Write( bOverwrite, 
                              strReadLengthsFile,
                              strReadFastbFile,
                              strReadQualbFile,
                              strReadTrimsFile,
                              strReadFastnFile );
}
                                  

void
ReadDataManager::SetPrefetchStrategy( PrefetchStrategy *pStrategy ) 
{
    mpReadSequenceMgr->SetPrefetchStrategy( pStrategy );
    mpReadNameMgr->SetPrefetchStrategy( pStrategy );
}

#include "assembly/BinaryVecDataManagerTemplate.h"
template class BinaryVecDataManager<Bool>;
