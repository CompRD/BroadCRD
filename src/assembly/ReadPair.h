///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// A ReadPair is a pair of two reads and their expected separation and
// the expected standard deviation of that separation.

#ifndef ASSEMBLY_READPAIR
#define ASSEMBLY_READPAIR

#include "assembly/ReadToken.h"
#include "assembly/ReadPairToken.h"
#include "assembly/ReadPairManager.h"

#include <utility>
#include <iostream>

using std::ostream;

class ReadPair
{
 public:
  ReadPair() {}
  
  ReadPair( const ReadPairToken &token )
    : mToken( token ) {}

  ReadPair( ReadPairManager *pRPM, const int id )
    : mToken( pRPM, id ) {}

  operator ReadPairToken() const { return mToken; }

  bool operator< ( const ReadPair &other ) const
  { return ( mToken < other.mToken ); }
  
  bool operator== ( const ReadPair &other ) const
  { return ( mToken == other.mToken ); }
  
  bool operator> ( const ReadPair &other ) const
  { return ( mToken > other.mToken ); }
  
  bool operator!= ( const ReadPair &other ) const
  { return ! ( *this == other ); }
  
  // Returns true if read exists in some assembly.
  bool IsValid( ) const { return mToken.IsValid(); }
  
  int  GetId( ) const { return mToken.GetId(); }
  

  pair<ReadToken,ReadToken> GetReads() const;
  
  ReadToken GetRead1() const;
  ReadToken GetRead2() const;
  ReadToken GetOtherRead( const ReadToken &oneRead ) const;
  
  int GetExpectedSep()        const;
  int GetExpectedStdDev()     const;

  int GetExpectedInsertSize() const;

  void SetExpectedSep( const int sep );
  void SetExpectedStdDev( const int stdev );
  
 private:
  ReadPairManager * GetMgrPtr() const { return mToken.GetMgrPtr(); }

  ReadPairToken mToken;
};


inline
ReadToken
ReadPair::GetRead1() const
{
  return GetMgrPtr()->GetRead1( mToken.GetId() );
}

inline
ReadToken
ReadPair::GetRead2() const
{
  return GetMgrPtr()->GetRead2( mToken.GetId() );
}

inline
ReadToken
ReadPair::GetOtherRead( const ReadToken &oneRead ) const
{
  ReadToken read1 = this->GetRead1();
  ReadToken read2 = this->GetRead2();
  
  if ( read1 == oneRead )
    return read2;
  else if ( read2 == oneRead )
    return read1;
  else
    return ReadToken();
}
  
inline
int 
ReadPair::GetExpectedSep() const
{
  return GetMgrPtr()->GetExpectedSep( mToken.GetId() );
}

inline
int 
ReadPair::GetExpectedStdDev() const
{
  return GetMgrPtr()->GetExpectedStdDev( mToken.GetId() );
}


inline
void
ReadPair::SetExpectedSep( const int sep )
{
  GetMgrPtr()->SetExpectedSep( mToken.GetId(), sep );
}

inline
void 
ReadPair::SetExpectedStdDev( const int stdev )
{
  GetMgrPtr()->SetExpectedStdDev( mToken.GetId(), stdev );
}


inline
ostream & operator<< ( ostream & out, const ReadPair &aReadPair )
{
  return out << aReadPair.GetRead1()
             << "-[" << aReadPair.GetExpectedSep() << "~" << aReadPair.GetExpectedStdDev() << "]-"
             << aReadPair.GetRead2();
}

#endif
