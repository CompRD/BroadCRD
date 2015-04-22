///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// A ReadLocationInSuper is the position and orientation of a given
// Read in a given Super.

#ifndef ASSEMBLY_READ_LOCATION_IN_SUPER
#define ASSEMBLY_READ_LOCATION_IN_SUPER

#include "system/Assert.h"

#include "assembly/ReadLocationInContig.h"
#include "assembly/ContigLocation.h"

#include <functional>

class ReadLocationInSuper
  : public Location<SuperToken,ReadToken>
{
 public:
  ReadLocationInSuper()
    : Location<SuperToken,ReadToken>()
  {}

  ReadLocationInSuper( const Location<SuperToken,ReadToken> &token )
    : Location<SuperToken,ReadToken>( token )
  {}

  // Create with explicit Super, Read, Interval, and Orientation.
  ReadLocationInSuper( const SuperToken &theSuper, const ReadToken &theRead, 
                       const Interval theInterval, const Orientation theOrientation )
    : Location<SuperToken,ReadToken>( theSuper, theRead, theInterval, theOrientation )
  {}

  // Create with Super, Read, Interval, and Orientation calculated
  // from a given ContigLocation and ReadLocation.  (The Contig
  // referenced by the two Locations must be the same.)
  ReadLocationInSuper( const ContigLocation &theContigLoc,
                       const ReadLocation &theReadLoc );

  SuperToken  GetSuper( ) const { return m_container; }
  ReadToken   GetRead( )  const { return m_element; }

  ReadLocationInSuper GetReversed() const;

  ContigLocation GetImpliedContigLocation( const ReadLocation &theReadLoc ) const;
};

inline
ReadLocationInSuper::ReadLocationInSuper( const ContigLocation &theContigLoc,
                                          const ReadLocation &theReadLoc )
  : Location<SuperToken,ReadToken>( theContigLoc.GetSuper(), theReadLoc.GetRead(),
                                    Interval(), orient_FW )
{
  ForceAssert( theContigLoc.GetContig() == theReadLoc.GetContig() );
  
  if ( theContigLoc.GetOrientation() == orient_FW )
  { 
    m_interval = Interval( theContigLoc.Begin() + theReadLoc.Begin(),
                           theContigLoc.Begin() + theReadLoc.End() );
    m_orientation = theReadLoc.GetOrientation();
  }
  else
  {
    m_interval = Interval( theContigLoc.End() - theReadLoc.End(),
                           theContigLoc.End() - theReadLoc.Begin() );
    m_orientation = Flip( theReadLoc.GetOrientation() );
  }
}

#endif
