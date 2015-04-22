///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// A ContigLocation is the location and orientation of a given Contig
// within a given Super.

#ifndef ASSEMBLY_CONTIG_LOCATION
#define ASSEMBLY_CONTIG_LOCATION

#include "assembly/Location.h"

#include "assembly/SuperToken.h"
#include "assembly/ContigToken.h"

#include <functional>


using std::binary_function;


class ContigLocation
  : public Location<SuperToken,ContigToken>
{
 public:
  ContigLocation( ) 
    : Location<SuperToken,ContigToken>() 
  {}

  ContigLocation( const Location<SuperToken,ContigToken> &token )
    : Location<SuperToken,ContigToken>( token )
  {}

  ContigLocation( const SuperToken &theSuper, const ContigToken &theContig, 
                  const Interval theInterval, const Orientation theOrientation )
    : Location<SuperToken,ContigToken>( theSuper, theContig, theInterval, theOrientation )
  {}

  ContigLocation( const SuperToken &theSuper, const ContigToken &theContig, 
                  const int begin, const Orientation theOrientation );

  SuperToken  GetSuper( )  const { return m_container; }
  ContigToken GetContig( ) const { return m_element; }

  ContigLocation GetReversed() const;
  
  // Functors.

  struct OrderBySuperAndLeastGap
    : public binary_function<ContigLocation,ContigLocation,bool>
  {
    bool
    operator() ( const ContigLocation &lhs,
                 const ContigLocation &rhs );
  };
  
  struct OrderPtrsBySuper
    : public binary_function<const ContigLocation*,const ContigLocation*,bool>
  {
    bool
    operator() ( const ContigLocation *pLHS,
                 const ContigLocation *pRHS );
  };
  
  struct OrderPtrsBySuperAndInterval
    : public binary_function<const ContigLocation*,const ContigLocation*,bool>
  {
    bool
    operator() ( const ContigLocation *pLHS,
                 const ContigLocation *pRHS );
  };
  
  struct OrderPtrsByContig 
    : public binary_function<const ContigLocation*,const ContigLocation*,bool>
  {
    bool
    operator() ( const ContigLocation *pLHS,
                 const ContigLocation *pRHS );
  };
};


inline bool
ContigLocation::OrderBySuperAndLeastGap::
operator() ( const ContigLocation &lhs,
	     const ContigLocation &rhs )
{
  return ( lhs.GetSuper().GetId() < rhs.GetSuper().GetId() ||
	   ( lhs.GetSuper().GetId() == rhs.GetSuper().GetId() &&
	     lhs.Begin() - rhs.End() < rhs.Begin() - lhs.End()) );
}


inline bool
ContigLocation::OrderPtrsBySuper::
operator() ( const ContigLocation *pLHS,
	     const ContigLocation *pRHS )
{
  return ( pLHS->GetSuper().GetId() < pRHS->GetSuper().GetId() );
}


inline bool
ContigLocation::OrderPtrsBySuperAndInterval::
operator() ( const ContigLocation *pLHS,
	     const ContigLocation *pRHS )
{
  return ( *pLHS < *pRHS );
}


inline bool
ContigLocation::OrderPtrsByContig::
operator() ( const ContigLocation *pLHS,
	     const ContigLocation *pRHS )
{
  return ( pLHS->GetContig().GetId() < pRHS->GetContig().GetId() );
}


#endif
