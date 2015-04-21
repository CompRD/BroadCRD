///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSEMBLY_READTOKEN
#define ASSEMBLY_READTOKEN

#include "assembly/Token.h"

class ReadDataManager;
class Read;

typedef Token<Read,ReadDataManager> ReadToken;

#include <iostream>

using std::ostream;

inline
ostream & operator<< ( ostream & out, const ReadToken &aReadToken )
{
    return out << 'R' << aReadToken.GetId();
}

#endif
