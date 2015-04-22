///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSEMBLY_READ
#define ASSEMBLY_READ

#include "assembly/ReadLocationInContig.h"
#include "assembly/ReadPair.h"
#include "assembly/ReadToken.h"
#include "assembly/SequenceIterator.h"

#include "Basevector.h"
#include "CompressedSequence.h"
#include "Qualvector.h"

#include <utility>

class Read
{
  public:
    Read( ) { }

    Read( ReadDataManager *pReadDataMgr, const int id )
        : mToken( pReadDataMgr, id ) { }

    Read( const ReadToken &readToken )
        : mToken( readToken ) { }

    operator ReadToken( ) const { return mToken; }

    bool operator< ( const Read &other ) const
    { return ( mToken < other.mToken ); }

    bool operator== ( const Read &other ) const
    { return ( mToken == other.mToken ); }

    bool operator> ( const Read &other ) const
    { return ( mToken > other.mToken ); }

    bool operator!= ( const Read &other ) const
    { return ! ( *this == other ); }

    // Returns true if read exists in some assembly.
    bool IsValid( ) const { return mToken.IsValid(); }

    int  GetId( ) const { return mToken.GetId(); }

    String                     GetName( )           const;
    int                        GetLength( )         const;
    const basevector &         GetBases( )          const;
    const qualvector &         GetQuals( )          const;
    const CompressedSequence & GetUntrimmedBases( ) const;
    pair<int,int>              GetTrims( )          const;
    int                        GetLeftTrim( )       const;
    int                        GetRightTrim( )      const;

    Bool IsRepetitive() const;

    // Methods to create SequenceIterators for a Read.
    SequenceIterator<Read> 
    GetSequenceBegin( Orientation orient )
    { return GetBeginOfSequence( *this, orient ); }

    SequenceIterator<Read> 
    GetSequenceEnd( Orientation orient)
    { return GetEndOfSequence( *this, orient ); }


    // Returns the number of placements of this read.
    int GetNumLocations( ) const;

    // Fills out vecRLs with all the ReadLocations of this Read.
    void GetLocations( vec<ReadLocation> &vecRLs ) const;

    Read     GetPartner( ) const; // Returns an invalid read if this read is unpaired.
    ReadPair GetPair   ( ) const; // Returns an invalid pair if this read is unpaired.

  private:
    ReadDataManager * GetMgrPtr() const { return mToken.GetMgrPtr(); }

    ReadToken mToken;
};


#endif
