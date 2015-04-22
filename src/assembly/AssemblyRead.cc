///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "assembly/AssemblyRead.h"

#include "assembly/ReadDataManager.h"

String
Read::GetName() const
{
    return GetMgrPtr()->GetName( GetId() );
}

int
Read::GetLength() const
{
    return GetMgrPtr()->GetLength( GetId() );
}

const basevector & 
Read::GetBases() const
{
    return GetMgrPtr()->GetBases( GetId() );
}

const qualvector & 
Read::GetQuals() const
{
    return GetMgrPtr()->GetQuals( GetId() );
}

const CompressedSequence & 
Read::GetUntrimmedBases() const
{
    return GetMgrPtr()->GetUntrimmedBases( GetId() );
}

pair<int,int> 
Read::GetTrims() const
{
    return GetMgrPtr()->GetTrims( GetId() );
}

int 
Read::GetLeftTrim() const
{
    return GetMgrPtr()->GetTrims( GetId() ).first;
}

int 
Read::GetRightTrim() const
{
    return GetMgrPtr()->GetTrims( GetId() ).second;
}

Bool
Read::IsRepetitive() const
{
    return GetMgrPtr()->IsRepetitive( GetId() );
}

int
Read::GetNumLocations() const
{
    return GetMgrPtr()->GetNumLocations( GetId() );
}

void
Read::GetLocations( vec<ReadLocation> &vecRLs ) const
{
    GetMgrPtr()->GetLocations( GetId(), vecRLs );
}

ReadPair
Read::GetPair( ) const
{
    return GetMgrPtr()->GetPair( GetId() );
}

Read
Read::GetPartner( ) const
{
  ReadPair thisPair = this->GetPair();
  if ( thisPair.IsValid() )
    return thisPair.GetOtherRead( mToken );
  else
    return Read();
}
