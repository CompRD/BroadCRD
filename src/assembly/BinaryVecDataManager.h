///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSEMBLY_BINARYVECDATAMANAGER
#define ASSEMBLY_BINARYVECDATAMANAGER

#include <utility>

#include "String.h"
#include "Vec.h"

template <class T>
class BinaryVecDataManager 
{
  public:
    BinaryVecDataManager( )
        : mbDataLoaded( true ),
          mbModified( false )
    { }

    BinaryVecDataManager( const String &strVecFile )
        : mStrVecFile( strVecFile ),
          mbDataLoaded( false ),
          mbModified( false )
    {
    }

    void LoadAllData( );

    void LoadData( const vec<int> &vecIds )
    {
        LoadAllData();
    }

    void LoadData( const int id )
    {
        LoadAllData();
    }

    const T & GetData( const int id )
    {
        if ( ! mbDataLoaded )
            LoadData( id );

        return mVecData[ id ];
    }

    void SetData( const int id, const T &data )
    {
        mbModified = true;
        LoadData( id );
        ResizeToFit( id );
        mVecData[ id ] = data;
    }

    unsigned int Size()
    {
        LoadAllData();
        return mVecData.size();
    }

    void Write( const bool bOverwrite,
                const String &strVecFile );

    bool fileExists() { return !mStrVecFile.empty(); }

  private:
    String mStrVecFile;

    bool mbDataLoaded;
    bool mbModified;

    vec<T> mVecData;
    
    void ResizeToFit( const int id );
};
extern template class BinaryVecDataManager<int>;

#endif
