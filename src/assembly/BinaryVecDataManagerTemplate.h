///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file BinaryVecDataManagerTemplate.h
 * \author tsharpe
 * \date Jun 27, 2012
 *
 * \brief
 */
#ifndef ASSEMBLY_BINARYVECDATAMANAGERTEMPLATE_H_
#define ASSEMBLY_BINARYVECDATAMANAGERTEMPLATE_H_

#include "assembly/BinaryVecDataManager.h"
#include "feudal/BinaryStream.h"

template <class T>
void
BinaryVecDataManager<T>::LoadAllData()
{
    if ( ! mbDataLoaded )
    {
        BinaryReader::readFile( mStrVecFile, &mVecData );
        mbDataLoaded = true;
    }
}

template <class T>
void
BinaryVecDataManager<T>::ResizeToFit( const int id )
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
BinaryVecDataManager<T>::Write( const bool bOverwrite,
                                const String &strOutFile )
{
    if ( ! mbModified &&
         IsRegularFile(strOutFile) &&
         IsRegularFile(mStrVecFile) &&
         RealPath(strOutFile) == RealPath(mStrVecFile) )
        return;

    ForceAssert( bOverwrite || ! mbModified || ! IsRegularFile( strOutFile ) );

    this->LoadAllData();

    BinaryWriter::writeFile( strOutFile, mVecData );
}

#endif /* ASSEMBLY_BINARYVECDATAMANAGERTEMPLATE_H_ */
