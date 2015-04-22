///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// The class  NameManager  defines the corresponding manger object
// which can be used to manage a collection of names (for reads, flows, etc).


#ifndef ASSEMBLY_NAMEMANAGER
#define ASSEMBLY_NAMEMANAGER

#include "String.h"
#include "Vec.h"

#include "assembly/FeudalDataManager.h"
#include "assembly/AddBehavior.h"

#include <map>

class PrefetchStrategy;

class NameManager 
{
  public:
    NameManager();
    NameManager( const String &strNameFile );
    ~NameManager();

    String GetName( const int id ) const;
    void   SetName( const int id, const String& name );

    int  FindId( const String &name ) const;

    bool Verify( const int id, AddBehavior b = ASSERT_IF_DIFFERENT ) const;

    
    void Write( const bool bOverwrite, 
                const String &NamesFile );


    // Give name manager a new prefetch strategy.  The name
    // manager will delete its old strategy and take responsibility
    // for managing the new one.
    void SetPrefetchStrategy( PrefetchStrategy *pPrefetchStrategy );

  private:
    StringsManager * mpStringsMgr;

    mutable bool mbDirty;
    mutable map<String,int> mMapNameToId;

    bool IsDirty() const;
    void Clean() const;

    PrefetchStrategy *mpPrefetchStrategy;
};
  
#endif

