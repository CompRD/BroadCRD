///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef GET_HOMES_H
#define GET_HOMES_H

#include "Basevector.h"
#include "CoreTools.h"
#include "paths/Uniseq.h"

void GetHomes( const String run_dir, const int K2, const vecbasevector& unibases2, 
     snark& S, const Bool VERBOSE );

#endif
