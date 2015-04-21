///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSEMBLY_CONTIGSEQUENCEMANAGER
#define ASSEMBLY_CONTIGSEQUENCEMANAGER

#include "Basevector.h"
#include "Qualvector.h"
#include "String.h"

#include "assembly/FeudalDataManager.h"

class PrefetchStrategy;

class ContigSequenceManager 
{
  public:
    ContigSequenceManager( );

    ContigSequenceManager( const String &strContigFastbFile,
                           const String &strContigQualbFile );

    ~ContigSequenceManager( );

    void LoadBases( const int id ) const;
    void LoadQuals( const int id ) const;

    void LoadBases( const vec<int> vecIds ) const;
    void LoadQuals( const vec<int> vecIds ) const;

    void LoadAllBases() const;
    void LoadAllQuals() const;

    const basevector & GetBases( const int id ) const;
    const qualvector & GetQuals( const int id ) const;

    void SetBases( const int id, const basevector &bases );
    void SetQuals( const int id, const qualvector &quals );

    bool Verify( const int id ) const;

    int GetLength( const int id ) const;

    void Reverse( const int id );


    void Write( const bool bOverwrite,
                const int lastId,
                const String &strContigFastbFile,
                const String &strContigQualbFile );

    // Give sequence manager a new prefetch strategy.  The sequence
    // manager will delete its old strategy and take responsibility
    // for managing the new one.
    void SetPrefetchStrategy( PrefetchStrategy *pPrefetchStrategy );

  private:
    BasesManager *mpBasesMgr;
    QualsManager *mpQualsMgr;

    PrefetchStrategy *mpPrefetchStrategy;
};
  
#endif
