///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "system/Assert.h"

#include "assembly/ReadPair.h"
#include "assembly/ReadPairManager.h"
#include "assembly/AssemblyRead.h"

pair<ReadToken,ReadToken> 
ReadPair::GetReads() const
{
  return make_pair( GetMgrPtr()->GetRead1( mToken.GetId() ),
                    GetMgrPtr()->GetRead2( mToken.GetId() ) );
}

int
ReadPair::GetExpectedInsertSize() const
{
    Read read1( this->GetRead1() );
    Read read2( this->GetRead2() );

    return ( read1.GetLeftTrim() + read1.GetLength() +
             this->GetExpectedSep() +
             read2.GetLength() + read2.GetLeftTrim() );
}

