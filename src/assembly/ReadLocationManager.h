///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSEMBLY_READLOCATIONMANAGER
#define ASSEMBLY_READLOCATIONMANAGER

#include "Vec.h"

#include "assembly/ReadLocationInContig.h"

#include <map>

class ContigDataManager;
class ReadDataManager;

class ReadLocationManager
{
  public:
    ReadLocationManager( ContigDataManager * &rpContigDataMgr,
                         ReadDataManager * &rpReadDataMgr,
                         const String &strReadLocationFile = "" );

    void Add( const ReadLocation &theReadLoc );

    bool Remove( const ReadLocation &theReadLoc );

    void ClearContig( const int contigId );

    void Reverse( const ContigToken theContig );

    void GetByContig( const int contigId, vec<ReadLocation> &vecReadLocs );

    void GetByRead( const int readId, vec<ReadLocation> &vecReadLocs );

    int  GetNumReadLocationsOfRead( const int readId );

    int  GetNumReadLocationsInContig( const int contigId );

    void Write( const bool bOverwrite,
                const String &strReadLocationFile );

  private:
    bool mbLoaded;
    bool mbModified;

    String mStrReadLocationFile;

    ContigDataManager * &mrpContigDataMgr;
    ReadDataManager * &mrpReadDataMgr;

    void Load();

    vec<ReadLocation> mVecReadLocations;
    vec< vec<int> >   mIndicesByContig;
    vec<int>          mIndicesByRead;
    multimap<int,int> mMapMultiPlacedReads;
    
    bool mbByContigIsOutOfDate;
    bool mbByReadIsOutOfDate;

    void UpdateByContig();
    void UpdateByRead();
};

#endif
