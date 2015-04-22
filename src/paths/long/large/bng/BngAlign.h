///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef BNG_ALIGN_H
#define BNG_ALIGN_H

#include "CoreTools.h"
#include "graph/Digraph.h"
#include "math/Functions.h"

// Align a BNG map R to a DISCOVAR assembly line A.  We assume that position r_m on
// R maps to position a_m on A.  This uses a shortest path algorithm.  Diagnostic
// output goes to "out".

double BngAlign( const vec<int>& R, const vec<int>& A, const int r_m, const int a_m,
     ostringstream& out, vec<int>& p, int& rstart, int& rstop, int& ndirect, 
     int& rsum, int& asum, const int mode );

#endif
