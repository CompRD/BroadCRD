///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "assembly/AssemblyContig.h"

ReadLocation
ReadLocation::GetReversed() const
{
  Contig theContig( this->GetContig() );
  
  return ReadLocation( theContig, this->GetRead(),
                       Interval( theContig.GetLength() - m_interval.End(),
                                 theContig.GetLength() - m_interval.Begin() ),
                       Flip( m_orientation ) );
}
