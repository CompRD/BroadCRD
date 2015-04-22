///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file VecDataManagerTemplate.h
 * \author tsharpe
 * \date Jun 27, 2012
 *
 * \brief
 */
#ifndef ASSEMBLY_VECDATAMANAGERTEMPLATE_H_
#define ASSEMBLY_VECDATAMANAGERTEMPLATE_H_

#include "assembly/VecDataManager.h"

template <class T>
void
VecDataManager<T>::LoadAllData()
{
    if ( ! mbDataLoaded )
    {
        Assert( IsAsciiVec( mStrVecFile ) );

        READX( mStrVecFile, mVecData );
        mbDataLoaded = true;
    }
}

template <class T>
void
VecDataManager<T>::ResizeToFit( const int id )
{
    unsigned int targetSize = id + 1;

    if ( mVecData.size() < targetSize )
    {
        unsigned int targetCapacity = mVecData.capacity();

        while ( targetCapacity < targetSize )
        {
            targetCapacity = (targetCapacity+1) * 2;
        }

        mVecData.reserve( targetCapacity );
        mVecData.resize( targetSize );
    }
}

template <class T>
void
VecDataManager<T>::Write( const bool bOverwrite,
                          const String &strOutFile )
{
    if ( ! mbModified &&
         IsRegularFile(strOutFile) &&
         IsRegularFile(mStrVecFile) &&
         RealPath(strOutFile) == RealPath(mStrVecFile) )
        return;

    ForceAssert( bOverwrite || ! mbModified || ! IsRegularFile( strOutFile ) );

    this->LoadAllData();

    WRITE( strOutFile, mVecData );
}

#endif /* ASSEMBLY_VECDATAMANAGERTEMPLATE_H_ */
