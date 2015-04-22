///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSEMBLY_READDATAMANAGER
#define ASSEMBLY_READDATAMANAGER

#include "Basevector.h"
#include "CompressedSequence.h"
#include "Qualvector.h"
#include "String.h"

#include "assembly/AssemblyRead.h"
#include "assembly/AddBehavior.h"
#include "assembly/BinaryVecDataManager.h"
#include <utility>

extern template class BinaryVecDataManager<Bool>;

class IdManager;
class ReadLocationManager;
class NameManager;
class ReadPairManager;
class ReadSequenceManager;
class PrefetchStrategy;

class ReadDataManager 
{
  public:
    ReadDataManager( ReadLocationManager *pRLM, 
                     ReadPairManager *pRPM );

    ReadDataManager( ReadLocationManager *pRLM, 
                     ReadPairManager *pRPM,
                     NameManager *pRNM,
                     ReadSequenceManager *pRSM,
                     BinaryVecDataManager<Bool> *pRFM,
                     IdManager *pIM );

    ~ReadDataManager( );


    void SetLocationManager( ReadLocationManager *pRLM );


    int GetSize( ) const;


    Read NewRead( const String &name,
                  const basevector &bases,
                  const qualvector &quals,
                  const pair<int,int> trims,
                  const CompressedSequence &untrimmedBases,
		  AddBehavior behavior = ASSERT_IF_DIFFERENT );

    Read GetRead     ( const int id );
    Read GetReadNamed( const String &name );

    String                     GetName          ( const int id ) const;
    int                        GetLength        ( const int id ) const;
    const basevector &         GetBases         ( const int id ) const;
    const qualvector &         GetQuals         ( const int id ) const;
    const CompressedSequence & GetUntrimmedBases( const int id ) const;
    
    pair<int,int>              GetTrims         ( const int id ) const;

    ReadPair                   GetPair          ( const int id ) const;

    Bool IsRepetitive( const int id ) const;

    int  GetNumLocations( const int id ) const;
    void GetLocations( const int id, vec<ReadLocation> &vecRLs ) const;


    void Write( const bool bOverwrite,
                const String &strReadNamesFile,
                const String &strReadLengthsFile,
                const String &strReadFastbFile,
                const String &strReadQualbFile,
                const String &strReadTrimsFile,
                const String &strReadFastnFile );

    void SetPrefetchStrategy( PrefetchStrategy *pStrategy );

  private:
    IdManager * mpIdMgr;
    NameManager * mpReadNameMgr;
    ReadSequenceManager * mpReadSequenceMgr;
    BinaryVecDataManager<Bool> * mpRepetitiveFlagMgr;

    ReadLocationManager * mpReadLocationMgr;
    ReadPairManager * mpReadPairMgr;

    bool CheckId( const int id ) const;
};
  
#endif

