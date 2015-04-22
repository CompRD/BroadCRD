///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "assembly/ContigLocation.h"

#include "assembly/AssemblyContig.h"

ContigLocation::ContigLocation( const SuperToken &theSuper, const ContigToken &theContig, 
                                const int begin, 
                                const Orientation theOrientation )
  : Location<SuperToken,ContigToken>( theSuper, theContig, 
                                      Interval( begin, begin + Contig(theContig).GetLength() ),
                                      theOrientation )
{}

#include "assembly/Super.h"

ContigLocation
ContigLocation::GetReversed() const
{
  Super theSuper( this->GetSuper() );
  
  return ContigLocation( theSuper, this->GetContig(),
                         Interval( theSuper.GetLength() - m_interval.End(),
                                   theSuper.GetLength() - m_interval.Begin() ),
                         Flip( m_orientation ) );
}

