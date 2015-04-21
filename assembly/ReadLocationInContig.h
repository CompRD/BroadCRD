///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// A ReadLocation is the position and orientation of a given Read in a
// given Contig.

#ifndef ASSEMBLY_READ_LOCATION_IN_CONTIG
#define ASSEMBLY_READ_LOCATION_IN_CONTIG

#include "assembly/Location.h"

#include "assembly/ReadToken.h"
#include "assembly/ContigToken.h"

#include <functional>


using std::binary_function;


class ReadLocation
  : public Location<ContigToken,ReadToken>
{
 public:
  ReadLocation()
    : Location<ContigToken,ReadToken>()
  {}

  ReadLocation( const Location<ContigToken,ReadToken> &token )
    : Location<ContigToken,ReadToken>( token )
  {}

  ReadLocation( const ContigToken &theContig, const ReadToken &theRead, 
                const Interval theInterval, const Orientation theOrientation )
    : Location<ContigToken,ReadToken>( theContig, theRead, theInterval, theOrientation ) 
  {}

  ContigToken GetContig( ) const { return m_container; }
  ReadToken   GetRead( )   const { return m_element; }

  ReadLocation GetReversed() const;


  // Modify 
  bool FlipOrientation() {
    if (m_orientation == orient_FW) {
      m_orientation = orient_RC;
      return true;
    }
    if (m_orientation == orient_RC) {
      m_orientation = orient_FW;
      return true;
    }
    return false; // presumably this line is not executed
  }
  

  // Functors.
  
  struct OrderPtrsByContig 
    : public binary_function<const ReadLocation*,const ReadLocation*,bool>
  {
    inline bool
    operator() ( const ReadLocation *pLHS,
                 const ReadLocation *pRHS );
  };
  
  
  struct EqualPtrsByContig 
    : public binary_function<const ReadLocation*,const ReadLocation*,bool>
  {
    inline bool
    operator() ( const ReadLocation *pLHS,
                 const ReadLocation *pRHS );
  };
  
  
  struct OrderPtrsByContigAndInterval 
    : public binary_function<const ReadLocation*,const ReadLocation*,bool>
  {
    inline bool
    operator() ( const ReadLocation *pLHS,
                 const ReadLocation *pRHS );
  };
  
  
  struct OrderPtrsByRead 
    : public binary_function<const ReadLocation*,const ReadLocation*,bool>
  {
    inline bool
    operator() ( const ReadLocation *pLHS,
                 const ReadLocation *pRHS );
  };
  
  
  struct EqualPtrsByRead 
    : public binary_function<const ReadLocation*,const ReadLocation*,bool>
  {
    inline bool
    operator() ( const ReadLocation *pLHS,
                 const ReadLocation *pRHS );
  };
  
};



inline bool
ReadLocation::OrderPtrsByContig::
operator() ( const ReadLocation *pLHS, 
	     const ReadLocation *pRHS )
{
  return ( pLHS->GetContig() < pRHS->GetContig() );
}

inline bool
ReadLocation::EqualPtrsByContig::
operator() ( const ReadLocation *pLHS,
	     const ReadLocation *pRHS )
{
  return ( pLHS->GetContig() == pRHS->GetContig() );
}

inline bool
ReadLocation::OrderPtrsByContigAndInterval::
operator() ( const ReadLocation *pLHS, 
	     const ReadLocation *pRHS )
{
  return ( pLHS->GetContig() < pRHS->GetContig() ||
	   ( pLHS->GetContig() == pRHS->GetContig() &&
	     pLHS->GetInterval() < pRHS->GetInterval() ) );
}

inline bool
ReadLocation::OrderPtrsByRead::
operator() ( const ReadLocation *pLHS, 
	     const ReadLocation *pRHS )
{
  return ( pLHS->GetRead() < pRHS->GetRead() );
}

inline bool
ReadLocation::EqualPtrsByRead::
operator() ( const ReadLocation *pLHS, 
	     const ReadLocation *pRHS )
{
  return ( pLHS->GetRead() == pRHS->GetRead() );
}



#endif
