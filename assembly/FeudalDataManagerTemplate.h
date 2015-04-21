///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// This file contains the implementation of the FeudalDataManager and
// is included into the locations that need to instantiate one.
#ifndef FEUDAL_DATA_MANAGER_TEMPLATE
#define FEUDAL_DATA_MANAGER_TEMPLATE

#include "assembly/FeudalDataManager.h"
#include "feudal/IncrementalWriter.h"

// TODO: Potentially dangerous truncation of ID
template <class MasterVecT, class SerfVecT>
void
FeudalDataManager<MasterVecT,SerfVecT>::LoadData( const vec<int> &vecIds )
{
    if ( ! mbInitialized )
        Initialize();

    vec<int> vecIdsToLoad;
    vecIdsToLoad.reserve( vecIds.size() );
    for ( unsigned int idIdx = 0; idIdx < vecIds.size(); ++idIdx )
    {
        int id = vecIds[ idIdx ];
        if ( ! mBvecDataLoaded[ id ] )
        {
            vecIdsToLoad.push_back( id );
            mBvecDataLoaded[ id ] = True;
        }
    }
    if ( vecIdsToLoad.size() )
        mVecData.SparseRead( mStrFeudalFile, vecIdsToLoad, 0 );
}

// TODO: Potentially dangerous truncation of ID
template <class MasterVecT, class SerfVecT>
void
FeudalDataManager<MasterVecT,SerfVecT>::SetData( const int id, const SerfVecT &data )
{
    if ( ! mbInitialized )
        Initialize();

    mbModified = true;
    ResizeToFit( id, data );

    // To trick the feudal vector into storing data internally, we push back data...
    mVecData.push_back( data );
    // swap it into its correct location...
    mVecData.SwapElements( id, mVecData.size() - 1 );
    // then remove the element that used to be at that index and is now the last element.
    mVecData.resize( mVecData.size() - 1 );

    mBvecDataLoaded[ id ] = True;
}

// TODO: Potentially dangerous truncation of ID
template <class MasterVecT, class SerfVecT>
void
FeudalDataManager<MasterVecT,SerfVecT>::ResizeToFit( unsigned long id,
                                                          const SerfVecT &data )
{
    size_t targetSize = id + 1;
    size_t targetCapacity = mVecData.capacity();

    const int startingMultiple = 10;

    // We use <= here to avoid always resizing in the push_back in SetData.
    if ( targetCapacity <= targetSize )
    {
        if ( targetCapacity == 0 )
            targetCapacity = startingMultiple;

        while ( targetCapacity <= targetSize )
            targetCapacity = targetCapacity * 2;
    }

    mVecData.reserve( targetCapacity );

    if ( static_cast<size_t>(mVecData.size()) < targetSize )
        mVecData.resize( targetSize );

    if ( mBvecDataLoaded.size() < targetSize )
        mBvecDataLoaded.resize( targetSize, False );
}

// TODO: Potentially dangerous truncation of ID
template <class MasterVecT, class SerfVecT>
void
FeudalDataManager<MasterVecT,SerfVecT>::Write( const bool bOverwrite,
                                               const String &strOutFile,
                                               const int lastId )
{
    if ( ! mbModified &&
         IsRegularFile(strOutFile) &&
         IsRegularFile(mStrFeudalFile) &&
         RealPath(strOutFile) == RealPath(mStrFeudalFile) )
      return;

    ForceAssert( bOverwrite || ! IsRegularFile( strOutFile ) );

    if ( ! mbModified && IsRegularFile( mStrFeudalFile ) ) {
      Cp( mStrFeudalFile, strOutFile );
      return;
    }

    if ( ! mbInitialized )
        Initialize();

    // Do I not need to worry about data I haven't read in yet?
    if ( mStrFeudalFile.empty() )
    {
        // All data is new, so just save existing vec.
        if ( lastId < 0 )
          mVecData.WriteAll( strOutFile );
        else
        {
          Remove( strOutFile );
          mVecData.WriteRange( strOutFile, 0, lastId+1 );
        }
    }
    else
    {
        // Some data may not have been loaded, so load in file in
        // pieces, saving to new location.
        MasterVecT tempVec;

        String outFile(strOutFile+".tmp");
        IncrementalWriter<SerfVecT> writer(outFile.c_str(), mBvecDataLoaded.size());

        int numPreexistingSerfs = MastervecFileObjectCount( mStrFeudalFile );

        if ( lastId >= 0 )
          numPreexistingSerfs = min( lastId + 1, numPreexistingSerfs );

        int chunkSize = max( 1, ( numPreexistingSerfs + 9 ) / 10 );

        int writeStart = 0, writeEnd = 0;
        for ( ; writeStart < numPreexistingSerfs; writeStart = writeEnd )
        {
          bool writeLocal = mBvecDataLoaded[ writeStart ];
          while ( writeEnd < numPreexistingSerfs &&
                  ( writeEnd - writeStart ) < chunkSize &&
                  mBvecDataLoaded[ writeEnd ] == writeLocal )
            ++writeEnd;

          if ( writeLocal )
              writer.add(mVecData.begin()+writeStart,mVecData.begin()+writeEnd);
          else
          {
            tempVec.ReadRange( mStrFeudalFile, writeStart, writeEnd );
            writer.add(tempVec.begin(),tempVec.end());
            tempVec.clear();
          }
        }

        // Write out any contigs we've added.
        int numExistingSerfs = mBvecDataLoaded.size();
        if ( lastId >= 0 )
          numExistingSerfs = min( lastId + 1, numExistingSerfs );

        if ( writeStart < numExistingSerfs )
            writer.add(mVecData.begin()+writeStart,mVecData.begin()+numExistingSerfs);

        writer.close();

        if ( IsRegularFile( strOutFile ) )
            Remove( strOutFile );

        Rename( strOutFile + ".tmp", strOutFile );
    }
}

#endif // FEUDAL_DATA_MANAGER_TEMPLATE
