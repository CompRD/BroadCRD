///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSEMBLY_CONTIGTOKEN
#define ASSEMBLY_CONTIGTOKEN

#include "assembly/Token.h"

class ContigDataManager;
class Contig;

typedef Token<Contig,ContigDataManager> ContigToken;

#include <iostream>

using std::ostream;

inline
ostream & operator<< ( ostream & out, const ContigToken &aContigToken )
{
    return out << 'C' << aContigToken.GetId();
}

#endif
