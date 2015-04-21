///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSEMBLY_TOKEN
#define ASSEMBLY_TOKEN

template <class ManagedObjectType, class MgrType>
class Token
{
  public:
    Token()
        : mpMgr( 0 ), mId( -1 ) { }

    Token( MgrType *pMgr, const int id )
        : mpMgr( pMgr ), mId( id ) { }

    typedef ManagedObjectType ObjectType;

    bool 
    operator== ( const Token &other ) const
    { return ( mpMgr == other.mpMgr && mId == other.mId ); }

    bool 
    operator!= ( const Token &other ) const
    { return ! ( *this == other ); }

    bool
    operator< ( const Token &other ) const
    { return ( mpMgr < other.mpMgr || mpMgr == other.mpMgr && mId < other.mId ); }

    bool
    operator> ( const Token &other ) const
    { return ( mpMgr > other.mpMgr || mpMgr == other.mpMgr && mId > other.mId ); }

    bool
    IsValid( ) const { return mpMgr != 0; }

    MgrType *  
    GetMgrPtr( ) const { return mpMgr; }

    int
    GetId( ) const { return mId; }

  private:
    MgrType *mpMgr;
    int mId;
};

#include <sys/types.h>

#if __GNUC__ > 2
namespace __gnu_cxx {
#endif

template <class T> struct hash;

template <class ManagedObjectType, class MgrType>
struct hash< Token<ManagedObjectType, MgrType> > {
    size_t 
    operator() ( const Token<ManagedObjectType, MgrType> &theToken ) const 
    { return theToken.GetId(); }
};

#if __GNUC__ > 2
}
#endif

#endif
