///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSEMBLY_FEUDALDATAMANAGER
#define ASSEMBLY_FEUDALDATAMANAGER

#include "Basevector.h"
#include "CompressedSequence.h"
#include "Qualvector.h"
#include "String.h"
#include "VecString.h"

// TODO: Potentially dangerous truncation of ID
// need to lose the vec<int>'s
template <class MasterVecT, class SerfVecT>
class FeudalDataManager
{
  public:
    FeudalDataManager( )
        : mbInitialized( true ),
          mbModified( false )
    { }

    FeudalDataManager( const String &strFeudalFile )
        : mStrFeudalFile( strFeudalFile ),
          mbInitialized( false ),
          mbModified( false )
    { }


    // Loads the entries in vecIds that have not already been loaded.
    // The vector vecIds must be sorted.
    void LoadData( const vec<int> &vecIds );

    void LoadData( const int id )
    {
        LoadData( vec<int>( 1, id ) );
    }

    void LoadAllData( )
    {
        int numIds = ( mStrFeudalFile.empty() ? 0 : MastervecFileObjectCount( mStrFeudalFile ) );

        vec<int> vecAllIds( numIds );
        for ( int id = 0; id < numIds; ++id )
            vecAllIds[ id ] = id;

        LoadData( vecAllIds );
    }

    const SerfVecT & GetData( const int id )
    {
        if ( ! IsLoaded( id ) )
            LoadData( id );

        return mVecData[ id ];
    }

    void SetData( const int id, const SerfVecT &data );

    void Write( const bool bOverwrite,
                const String &strFeudalFile,
                const int lastId = -1 );

    Bool IsLoaded( const int id )
    {
        if ( ! mbInitialized )
            Initialize();

        return mBvecDataLoaded[ id ];
    }

    // TODO: Potentially dangerous truncation of quals size
    int Size()
    {
        if ( ! mbInitialized )
            Initialize();

        return mVecData.size();
    }

  private:
    String mStrFeudalFile;

    bool mbInitialized;

    void Initialize( )
    {
        if ( mStrFeudalFile.size() )
        {
            int size = MastervecFileObjectCount( mStrFeudalFile );
            mVecData.reserve( size );
            mBvecDataLoaded.resize( size, False );
        }

        mbInitialized = true;
    }

    vec<Bool> mBvecDataLoaded;
    bool mbModified;

    MasterVecT mVecData;

    void ResizeToFit( unsigned long id, const SerfVecT &data );
};

typedef FeudalDataManager<vecbasevector,basevector> BasesManager;
typedef FeudalDataManager<vecqvec,qvec> QualsManager;
typedef FeudalDataManager<veccompseq,CompressedSequence> CompressedSequenceManager;
typedef FeudalDataManager<vecString,String> StringsManager;
extern template class FeudalDataManager<vecbasevector,basevector>;
extern template class FeudalDataManager<vecqvec,qvec>;
extern template class FeudalDataManager<veccompseq,CompressedSequence>;
extern template class FeudalDataManager<vecString,String>;

#endif
