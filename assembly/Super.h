///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// A Super may be thought of as a collection of ContigLocations.

#ifndef ASSEMBLY_SUPER
#define ASSEMBLY_SUPER

#include "Vec.h"
#include "assembly/ContigLocation.h"
#include "assembly/AssemblyLink.h"
#include "assembly/SuperToken.h"

class Super
{
  public:
    Super() { }

    Super( SuperDataManager *pSuperDataMgr, const int id )
        : mToken( pSuperDataMgr, id ) { }

    Super( const SuperToken &superToken )
        : mToken( superToken ) { }

    operator SuperToken() const { return mToken; }

    bool operator< ( const Super &other ) const
    { return ( mToken < other.mToken ); }

    bool operator== ( const Super &other ) const
    { return ( mToken == other.mToken ); }

    bool operator> ( const Super &other ) const
    { return ( mToken > other.mToken ); }

    bool operator!= ( const Super &other ) const
    { return ! ( *this == other ); }

    // Returns true if super exists in some assembly.
    bool IsValid( ) const { return mToken.IsValid(); }

    int  GetId( ) const { return mToken.GetId(); }

    // Returns distance from the lowest begin coordinate to highest
    // end coordinate of contained contigs.
    int  GetLength() const;

    // Returns sum of the lengths of the contained contigs.
    longlong  GetSumOfContigLengths() const;

    // Returns the number of ContigLocations in this Super.
    int  GetNumContigLocations( ) const;

    // Fills our vecCLs with all the ContigLocations in this Super.
    void GetContigLocations( vec<ContigLocation> &vecCLs ) const;

    // Returns the number of ReadLocations in this Super.
    int  GetNumReadLocations( ) const;

    // Fills out vecLinks will all the Links involving all the Reads
    // in all the Contigs in this Super.
    void GetLinks( vec<Link> &vecLinks ) const;

    void Print( ostream &out ) const;

  private:
    SuperDataManager * GetMgrPtr() const { return mToken.GetMgrPtr(); }

    SuperToken mToken;
};


#endif
