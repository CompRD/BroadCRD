///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Jumpster.  For each unibase that dead ends on the right, find its 40-mer
// that is back K from the end, and walk right through jump 40-mers for up to 1000
// bases.  At each step consider the 40+5-mers at that point; these must vote
// essentially unanimously.

#ifndef JUMPSTER_H
#define JUMPSTER_H

#include "Basevector.h"
#include "CoreTools.h"
#include "IntPairVec.h"
#include "Qualvector.h"

void Jumpster( const int K, const vecbasevector& unibases, 
     const vec< vec<int> >& nexts, const vecbasevector& jumps, 
     const vecqualvector& quals, vec<basevector>& bridges, const Bool verbose, 
     const vecbasevector& genome2, const int LG, 
     const VecIntPairVec& Glocs );

#endif
