///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSEMBLY_READPAIRTOKEN
#define ASSEMBLY_READPAIRTOKEN

#include "assembly/Token.h"

class ReadPairManager;
class ReadPair;

typedef Token<ReadPair,ReadPairManager> ReadPairToken;

#endif
