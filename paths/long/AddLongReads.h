///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#ifndef PATHS_LONG_ADDLONGREADS_H
#define PATHS_LONG_ADDLONGREADS_H
#include "Basevector.h"
#include "paths/long/SupportedHyperBasevector.h"

void UnrollWithPacbioReads( const vecbasevector& longreads, SupportedHyperBasevector * shb, 
        bool useOldLRPMethod, bool CorrectPacbioConsensus = true, int verbosity = 1);

#endif
