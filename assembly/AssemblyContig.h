///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// A Contig may be thought of as a collection of ReadLocations, but it
// is also a sequence.

#ifndef ASSEMBLY_CONTIG
#define ASSEMBLY_CONTIG

#include "Basevector.h"
#include "Qualvector.h"
#include "Vec.h"

#include "assembly/ContigLocation.h"
#include "assembly/ContigToken.h"
#include "assembly/AssemblyLink.h"
#include "assembly/ReadLocationInContig.h"
#include "assembly/SequenceIterator.h"

class Contig
{
  public:
    Contig( ) {}

    Contig( ContigDataManager *pContigDataMgr, const int id )
        : mToken( pContigDataMgr, id ) { }

    Contig( const ContigToken &contigToken )
        : mToken( contigToken ) { }

    operator ContigToken() const { return mToken; }

    bool operator< ( const Contig &other ) const
    { return ( mToken < other.mToken ); }

    bool operator== ( const Contig &other ) const
    { return ( mToken == other.mToken ); }

    bool operator> ( const Contig &other ) const
    { return ( mToken > other.mToken ); }

    bool operator!= ( const Contig &other ) const
    { return ! ( *this == other ); }

    // Returns true if contig exists in some assembly.
    bool IsValid( ) const { return mToken.IsValid(); }

    int  GetId( ) const { return mToken.GetId(); }

    // Returns length of consensus sequence.
    int  GetLength( ) const;
    
    // Returns the number of ReadLocations in this Contig.
    int  GetNumReadLocations( ) const;

    const basevector & GetBases( ) const;
    const qualvector & GetQuals( ) const;


    // Methods to create SequenceIterators for a Contig.
    SequenceIterator<Contig> 
    GetSequenceBegin( Orientation orient ) const
    { return GetBeginOfSequence( *this, orient ); }

    SequenceIterator<Contig> 
    GetSequenceEnd( Orientation orient ) const
    { return GetEndOfSequence( *this, orient ); }


    // Fills out vecRLs with all the ReadLocations in this Contig.
    void GetReadLocations( vec<ReadLocation> &vecRLs )   const;

    // Fills out vecCLs with all the ContigLocations of this Contig.
    void GetSelfLocations( vec<ContigLocation> &vecCLs ) const;

    // Fills out vecLinks will all the Links involving all the Reads
    // in this Contig.
    void GetLinks( vec<Link> &vecLinks ) const;

    // Fills out vecLinks will all the Links involving this Contig and
    // some other Contig.
    void GetExternalLinks( vec<Link> &vecLinks ) const;

  private:
    // Fills out vecLinks will all the Links involving all the Reads
    // in this Contig.
    void GetLinks( vec<Link> &vecLinks, bool externalOnly ) const;

    ContigDataManager * GetMgrPtr() const { return mToken.GetMgrPtr(); }
    
    ContigToken mToken;
};

#endif
