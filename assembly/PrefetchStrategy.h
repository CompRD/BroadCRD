///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSEMBLY_PREFETCHSTRATEGY
#define ASSEMBLY_PREFETCHSTRATEGY

#include "Vec.h"


// Abstract interface.

class PrefetchStrategy
{
  public:
    PrefetchStrategy( ) { }
    virtual ~PrefetchStrategy( ) { }

    virtual void GetIdsToLoad( const int seedId, vec<int> &vecIdsToLoad ) = 0;
  
    virtual bool GetAll() = 0;

    virtual PrefetchStrategy * Clone() const = 0;
};



// Returns only the id it's given.

class NullPrefetchStrategy : public PrefetchStrategy
{
  public:
    NullPrefetchStrategy( ) 
        : PrefetchStrategy( ) 
    { }
    
    virtual ~NullPrefetchStrategy( ) { }

    virtual void GetIdsToLoad( const int seedId, vec<int> &vecIdsToLoad )
    { 
        vecIdsToLoad = vec<int>( 1, seedId ); 
    }
  
    virtual bool GetAll() { return false; }

    virtual PrefetchStrategy * Clone() const { return new NullPrefetchStrategy( *this ); }
};



// Returns only the id it's given, but returns true for GetAll().

class PrefetchAllStrategy : public PrefetchStrategy
{
  public:
    PrefetchAllStrategy()
      : PrefetchStrategy( )
    { }

    virtual ~PrefetchAllStrategy( ) { }
     
    virtual void GetIdsToLoad( const int seedId, vec<int> &vecIdsToLoad )
    {
      vecIdsToLoad = vec<int>( 1, seedId );
    }

    virtual bool GetAll() { return true; }

    virtual PrefetchStrategy * Clone() const { return new PrefetchAllStrategy( *this ); }
};


class ReadLocationManager; 



// Returns the ids of all the reads in the same contig(s) as the read
// given, as well as their partners.

class PrefetchReadsByContigStrategy : public PrefetchStrategy
{
  public:
    PrefetchReadsByContigStrategy( ReadLocationManager *pReadLocationMgr,
                                   unsigned int neighborhood = 50 ) 
        : PrefetchStrategy( ),
          mpReadLocationMgr( pReadLocationMgr ),
          mNeighborhood( neighborhood )
    { }
    
    virtual ~PrefetchReadsByContigStrategy( ) { }

    virtual void GetIdsToLoad( const int seedId, vec<int> &vecIdsToLoad );

    virtual bool GetAll() { return false; }

    virtual PrefetchStrategy * Clone() const { return new PrefetchReadsByContigStrategy( *this ); }

  private:
    ReadLocationManager *mpReadLocationMgr;
    unsigned int mNeighborhood;
};



// Returns the ids of all the reads in the same super(s) as the read
// given, as well as their partners.

class PrefetchReadsBySuperStrategy : public PrefetchStrategy
{
  public:
    PrefetchReadsBySuperStrategy( ReadLocationManager *pReadLocationMgr,
                                  unsigned int neighborhood = 5 ) 
        : PrefetchStrategy( ),
          mpReadLocationMgr( pReadLocationMgr ),
          mNeighborhood( neighborhood )
    { }
    
    virtual ~PrefetchReadsBySuperStrategy( ) { }

    virtual void GetIdsToLoad( const int seedId, vec<int> &vecIdsToLoad );

    virtual bool GetAll() { return false; }

    virtual PrefetchStrategy * Clone() const { return new PrefetchReadsBySuperStrategy( *this ); }

  private:
    ReadLocationManager *mpReadLocationMgr;
    unsigned int mNeighborhood;
};



class ContigLocationManager; 


// Returns the ids of all the contigs in the same super(s) as the contig
// given.

class PrefetchContigsBySuperStrategy : public PrefetchStrategy
{
  public:
    PrefetchContigsBySuperStrategy( ContigLocationManager *pContigLocationMgr ) 
        : PrefetchStrategy( ),
          mpContigLocationMgr( pContigLocationMgr )
    { }
    
    virtual ~PrefetchContigsBySuperStrategy( ) { }

    virtual void GetIdsToLoad( const int seedId, vec<int> &vecIdsToLoad );

    virtual bool GetAll() { return false; }

    virtual PrefetchStrategy * Clone() const { return new PrefetchContigsBySuperStrategy( *this ); }

  private:
    ContigLocationManager *mpContigLocationMgr;
};

#endif
