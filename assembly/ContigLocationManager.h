///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSEMBLY_CONTIGLOCATIONMANAGER
#define ASSEMBLY_CONTIGLOCATIONMANAGER

#include "Vec.h"

#include "assembly/ContigLocation.h"

class ContigLocationManager
{
  public:
    ContigLocationManager( SuperDataManager * &rpSuperDataMgr,
                           ContigDataManager * &rpContigDataMgr );
    
    void Add( const ContigLocation &theContigLoc );

    bool Remove( const ContigLocation &theContigLoc );

    void ClearSuper( const int superId );

    void Reverse( const SuperToken theSuper );

    void Shift( const int superId, const int shiftAmount );

    void Normalize( const int superId );

  
    int  GetNumContigLocationsInSuper( const int superId );

    void GetBySuper( const int superId, vec<ContigLocation> &vecContigLocs );

    void GetByContig( const int contigId, vec<ContigLocation> &vecContigLocs );

    int  GetSuperLength( const int superId );

    longlong GetSumOfContigLengths( const int superId );

    bool WasModified( ) const { return mbModified; }

  private:
    SuperDataManager * &mrpSuperDataMgr;
    ContigDataManager * &mrpContigDataMgr;

    vec<ContigLocation> mVecContigLocations;

    vec< vec<int> > mIndicesBySuper;
    vec< vec<int> > mIndicesByContig;
    
    vec<int> mVecSuperBegin;
    vec<int> mVecSuperEnd;

    vec<longlong> mSumContigLengths;

    bool mbDirty;
    bool mbModified;

    unsigned int mNumEmpty;

    void Clean();
    
};

#endif
