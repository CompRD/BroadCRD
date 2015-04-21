///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef BNG_ALIGN_CORE_H
#define BNG_ALIGN_CORE_H

#include "CoreTools.h"

void BngAlignCore( const vec<vec<double>>& X,
     const vec< quad<int,int,int,int> >& zmatches,
     const int64_t low, const int64_t high,
     const vec<int>& S,
     const Bool ALIGN_LOGGING, const Bool SHOW_ALL_SEEDS, const Bool ALIGN_DETAILS,
     const Bool ONE_SEED, const Bool RFILTER, int& PRINT_FAIL,
     int& align_calls,
     vec< triple< int, int, vec<int> > >& aligns, vec<double>& scores,
     double& ac1, double& ac2, double& ac3, double& ac4, double& ac5,
     const int min_direct, const double max_score, const double max_total_err,
     const int s1, const int i1 );

#endif
