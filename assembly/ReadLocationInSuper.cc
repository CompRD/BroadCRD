///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "assembly/ReadLocationInSuper.h"

#include "assembly/AssemblyContig.h"
#include "assembly/Super.h"

ContigLocation
ReadLocationInSuper::GetImpliedContigLocation( const ReadLocation &theReadLoc ) const
{
  Contig theContig( theReadLoc.GetContig() );
  
  int contigBegin;
  Orientation contigOrient;
  
  if ( theReadLoc.GetOrientation() == m_orientation )
  {
    contigBegin = this->Begin() - theReadLoc.Begin();
    contigOrient = orient_FW;
  }
  else
  {
    contigBegin = this->Begin() - ( theContig.GetLength() - theReadLoc.End() );
    contigOrient = orient_RC;
  }
  
  int contigEnd = contigBegin + theContig.GetLength();
  
  return ContigLocation( this->GetSuper(), theContig, 
                         Interval( contigBegin, contigEnd ), contigOrient );
}

ReadLocationInSuper
ReadLocationInSuper::GetReversed() const
{
  Super theSuper( this->GetSuper() );
  
  return ReadLocationInSuper( theSuper, this->GetRead(),
                              Interval( theSuper.GetLength() - m_interval.End(),
                                        theSuper.GetLength() - m_interval.Begin() ),
                              Flip( m_orientation ) );
}
