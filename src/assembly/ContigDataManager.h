///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSEMBLY_CONTIGDATAMANAGER
#define ASSEMBLY_CONTIGDATAMANAGER

#include "Basevector.h"
#include "Qualvector.h"
#include "String.h"
#include "Vec.h"

#include "assembly/AssemblyContig.h"

class IdManager;
class ContigSequenceManager;
class ContigLocationManager;
class ReadLocationManager;
class PrefetchStrategy;

class ContigDataManager 
{
  public:
    ContigDataManager( ContigLocationManager *pCLM, 
                       ReadLocationManager *pRLM );

    ContigDataManager( ContigLocationManager *pCLM, 
                       ReadLocationManager *pRLM,
                       ContigSequenceManager *pCSM, 
                       IdManager *pIM );

    ~ContigDataManager( );


    void SetLocationManager( ContigLocationManager *pCLM );


    int GetSize( ) const;


    Contig NewContig( const basevector &bases,
                      const qualvector &quals );

    Contig GetContig( const int id );


    int GetLength( const int id ) const;

    int GetNumReadLocations( const int id ) const;

    const basevector & GetBases( const int id ) const;
    const qualvector & GetQuals( const int id ) const;

    void SetBases( const int id, const basevector &newBases );
    void SetQuals( const int id, const qualvector &newQuals );

    void GetReadLocations  ( const int id, vec<ReadLocation> &vecRLs )   const;
    void GetContigLocations( const int id, vec<ContigLocation> &vecCLs ) const;


    void Reverse( const int id );

    int PurgeEmptyContigs();

    void Write( const bool bOverwrite,
                const String &strContigFastbFile,
                const String &strContigQualbFile );

    void SetPrefetchStrategy( PrefetchStrategy *pStrategy );

  private:
    IdManager             * mpIdMgr;
    ContigSequenceManager * mpContigSequenceMgr;

    ContigLocationManager * mpContigLocationMgr;
    ReadLocationManager   * mpReadLocationMgr;

    bool CheckId( const int id ) const;
};
  
#endif

