///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSEMBLY_READSEQUENCEMANAGER
#define ASSEMBLY_READSEQUENCEMANAGER

#include "String.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "CompressedSequence.h"

#include "assembly/FeudalDataManager.h"
#include "assembly/VecDataManager.h"
#include "assembly/BinaryVecDataManager.h"

typedef BinaryVecDataManager<int> BinaryVecIntManager;

typedef VecDataManager< pair<int,int> > VecTrimsManager;
template <> void VecDataManager< pair<int,int> >::LoadAllData();
template <> void VecDataManager< pair<int,int> >::Write( const bool, const String & );
extern template class VecDataManager< std::pair<int,int> >;

class PrefetchStrategy;

class ReadSequenceManager 
{
  public:
    ReadSequenceManager( );

    ReadSequenceManager( const String &strReadLengthsFile,
                         const String &strReadFastbFile,
                         const String &strReadQualbFile,
                         const String &strReadTrimsFile,
                         const String &strReadFastnFile );

    ~ReadSequenceManager( );

    int                        GetLength        ( const int id ) const;
    const basevector &         GetBases         ( const int id ) const;
    const qualvector &         GetQuals         ( const int id ) const; 
    pair<int,int>              GetTrims         ( const int id ) const;
    const CompressedSequence & GetUntrimmedBases( const int id ) const;

    void SetBases         ( const int id, const basevector &bases );
    void SetQuals         ( const int id, const qualvector &quals );
    void SetTrims         ( const int id, const pair<int,int> trims );
    void SetUntrimmedBases( const int id, const CompressedSequence &untrimmedBases );

    bool Verify( const int id ) const;


    void Write( const bool bOverwrite, 
                const String &strReadLengthsFile,
                const String &strReadFastbFile,
                const String &strReadQualbFile,
                const String &strReadTrimsFile,
                const String &strReadFastnFile );

    // Give sequence manager a new prefetch strategy.  The sequence
    // manager will delete its old strategy and take responsibility
    // for managing the new one.
    void SetPrefetchStrategy( PrefetchStrategy *pPrefetchStrategy );

  private:
    BinaryVecIntManager *mpLengthsMgr;
    BasesManager  *mpBasesMgr;
    QualsManager  *mpQualsMgr;

    VecTrimsManager *mpTrimsMgr;

    CompressedSequenceManager *mpUntrimmedBasesMgr;

    PrefetchStrategy *mpPrefetchStrategy;
};
  
#endif

