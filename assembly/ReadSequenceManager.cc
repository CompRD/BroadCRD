///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "assembly/ReadSequenceManager.h"

#include "assembly/PrefetchStrategy.h"

ReadSequenceManager::ReadSequenceManager( )
    : mpLengthsMgr( new BinaryVecIntManager ),
      mpBasesMgr( new BasesManager ),
      mpQualsMgr( new QualsManager ),
      mpTrimsMgr( new VecTrimsManager ),
      mpUntrimmedBasesMgr( new CompressedSequenceManager ),
      mpPrefetchStrategy( new NullPrefetchStrategy )
{
}

ReadSequenceManager::ReadSequenceManager( const String &strReadLengthsFile,
                                          const String &strReadFastbFile,
                                          const String &strReadQualbFile,
                                          const String &strReadTrimsFile,
                                          const String &strReadFastnFile )
    : mpLengthsMgr( new BinaryVecIntManager( strReadLengthsFile ) ),
      mpBasesMgr( new BasesManager( strReadFastbFile ) ),
      mpQualsMgr( new QualsManager( strReadQualbFile ) ),
      mpTrimsMgr( new VecTrimsManager( strReadTrimsFile ) ),
      mpUntrimmedBasesMgr( new CompressedSequenceManager( strReadFastnFile ) ),
      mpPrefetchStrategy( new NullPrefetchStrategy )
{
}

ReadSequenceManager::~ReadSequenceManager( )
{
    delete mpLengthsMgr;
    delete mpBasesMgr;
    delete mpQualsMgr;
    delete mpTrimsMgr;
    delete mpUntrimmedBasesMgr;
    delete mpPrefetchStrategy;
}

int
ReadSequenceManager::GetLength( const int id ) const
{
    return mpLengthsMgr->GetData( id );
}

const basevector & 
ReadSequenceManager::GetBases( const int id ) const
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
ReadSequenceManager::GetQuals( const int id ) const
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

pair<int,int>
ReadSequenceManager::GetTrims( const int id ) const
{
    return mpTrimsMgr->GetData( id );
}

const CompressedSequence & 
ReadSequenceManager::GetUntrimmedBases( const int id ) const
{
    if ( ! mpUntrimmedBasesMgr->IsLoaded( id ) )
    {
        if ( mpPrefetchStrategy->GetAll() )
            mpUntrimmedBasesMgr->LoadAllData();
        else
        {
            vec<int> vecIds;
            mpPrefetchStrategy->GetIdsToLoad( id, vecIds );
            mpUntrimmedBasesMgr->LoadData( vecIds );
        }
    }

    return mpUntrimmedBasesMgr->GetData( id );
}


void
ReadSequenceManager::SetBases( const int id, const basevector &bases )
{
    mpBasesMgr->SetData( id, bases );
    mpLengthsMgr->SetData( id, bases.size() );
}

void
ReadSequenceManager::SetQuals( const int id, const qualvector &quals )
{
    mpQualsMgr->SetData( id, quals );
}

void
ReadSequenceManager::SetTrims( const int id, const pair<int,int> trims )
{
    mpTrimsMgr->SetData( id, trims );
}

void 
ReadSequenceManager::SetUntrimmedBases( const int id, const CompressedSequence &untrimmedBases )
{
    mpUntrimmedBasesMgr->SetData( id, untrimmedBases );
}


bool 
ReadSequenceManager::Verify( const int id ) const
{
    unsigned int basesLength = GetBases(id).size();

    ForceAssertEq( basesLength, GetQuals(id).size() );
    ForceAssertEq( basesLength, static_cast<unsigned>(GetLength(id)) );

    pair<int,int> trims = GetTrims(id);
    if ( trims.first + basesLength + trims.second !=
            static_cast<unsigned>(GetUntrimmedBases(id).size()) )
    {
        // This generates a lot of messages.  Oh, well.
        // cerr << "WARNING: left trim + trimmed length + right trim"
        //      << " (" << trims.first + basesLength + trims.second << ")"
        //      << " != untrimmed length"
        //      << " (" << GetUntrimmedBases(id).size() << ")"
        //      << " for read #" << id << endl;
    }

    return true;
}

void
ReadSequenceManager::Write( const bool bOverwrite, 
                            const String &strReadLengthsFile,
                            const String &strReadFastbFile,
                            const String &strReadQualbFile,
                            const String &strReadTrimsFile,
                            const String &strReadFastnFile )
{
    mpLengthsMgr->Write( bOverwrite, strReadLengthsFile );
    mpBasesMgr->Write  ( bOverwrite, strReadFastbFile );
    mpQualsMgr->Write  ( bOverwrite, strReadQualbFile );
    mpTrimsMgr->Write  ( bOverwrite, strReadTrimsFile );
    mpUntrimmedBasesMgr->Write( bOverwrite, strReadFastnFile );
}

void
ReadSequenceManager::SetPrefetchStrategy( PrefetchStrategy *pPrefetchStrategy )
{
    delete mpPrefetchStrategy;
    mpPrefetchStrategy = pPrefetchStrategy->Clone();
}

// specialize LoadAllData for trims.
template <>
void
VecDataManager< pair<int,int> >::LoadAllData( )
{
    if ( ! mbDataLoaded )
    {
        int nLines = LineCount( mStrVecFile );
        mVecData.reserve( nLines );

        Ifstream( trimsStream, mStrVecFile );

        for (int lineNum=0; lineNum<nLines; ++lineNum)
        {
            int leftAmount;
            int rightAmount;

            trimsStream >> leftAmount >> rightAmount;
            mVecData.push_back( make_pair( leftAmount, rightAmount ) );
        }

        trimsStream.close( );
        mbDataLoaded = true;
    }
}

// specialize Write for trims.
template <>
void
VecDataManager< pair<int,int> >::Write( const bool bOverwrite,
                                        const String &strOutFile )
{
    if ( ! mbModified &&
         IsRegularFile(strOutFile) &&
         IsRegularFile(mStrVecFile) &&
         RealPath(strOutFile) == RealPath(mStrVecFile) )
        return;

    ForceAssert( bOverwrite || ! IsRegularFile( strOutFile ) );

    this->LoadAllData();

    Ofstream( trimsStream, strOutFile );

    int nLines = mVecData.size();

    for (int lineNum=0; lineNum<nLines; ++lineNum)
    {
        pair<int,int> theTrims = mVecData[lineNum];

        trimsStream << theTrims.first << ' ' << theTrims.second << '\n';
    }

    trimsStream.close( );
}

#include "assembly/BinaryVecDataManagerTemplate.h"
template class BinaryVecDataManager<int>;

#include "assembly/VecDataManagerTemplate.h"
template class VecDataManager< pair<int,int> >;
