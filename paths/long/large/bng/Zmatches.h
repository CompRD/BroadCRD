///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ZMATCHES_H
#define ZMATCHES_H

#include "CoreTools.h"

// Build alignment seeds.

void Zmatches( const vec<vec<double>>& X, const vec<int>& S,
     vec< quad<int,int,int,int> >& zmatches, const Bool FILTER_DUPS,
     const int div, const int w, const int cutoff );

#endif
