///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSEMBLY_SUPERDATAMANAGER
#define ASSEMBLY_SUPERDATAMANAGER

#include "assembly/Super.h"
#include "assembly/ContigLocation.h"

class ContigLocationManager;
class IdManager;

class SuperDataManager 
{
  public:
    SuperDataManager( ContigLocationManager *pCLM );

    SuperDataManager( ContigLocationManager *pCLM,
                      IdManager *pIM );

    ~SuperDataManager( );


    int GetSize( ) const;


    Super NewSuper( );

    Super GetSuper( const int id );

    int  GetNumContigLocations( const int id ) const;

    void GetContigLocations( const int id, vec<ContigLocation> &vecCLs ) const;
    int  GetLength( const int id ) const;
    longlong GetSumOfContigLengths( const int id ) const;
    
    void Reverse( const int id );

  private:
    IdManager * mpIdMgr;

    ContigLocationManager * mpContigLocationMgr;
};
  
#endif

