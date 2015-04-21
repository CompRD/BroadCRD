///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSEMBLY_SUPERTOKEN
#define ASSEMBLY_SUPERTOKEN

#include "assembly/Token.h"

class SuperDataManager;
class Super;

typedef Token<Super,SuperDataManager> SuperToken;


#include <iostream>

using std::ostream;

inline
ostream & operator<< ( ostream & out, const SuperToken &aSuperToken )
{
    return out << 'S' << aSuperToken.GetId();
}
    


#endif
