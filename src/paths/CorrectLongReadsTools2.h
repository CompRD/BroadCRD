///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef CORRECT_LONG_READS_TOOLS2_H
#define CORRECT_LONG_READS_TOOLS2_H

#include "Basevector.h"
#include "MainTools.h"
#include "IntPairVec.h"
#include "graph/Digraph.h"
#include "paths/BigMapTools.h"
#include "paths/Useq.h"

void Assess( const vec< vec<int> >& reads, const int nr, const vecbasevector& genome,
     const vecbasevector& genome2, const vecbasevector& unibases, const int K, 
     const int LG, const VecIntPairVec& Glocs,
     const vec<int>& uperfect );

void CleanBubbles( const int K, const vecbasevector& unibases, 
     const vec< vec<int> >& nexts, vec< vec< vec<int> > >& upaths, 
     vec< vec<int> >& cores, vec< vec<int> >& ucores, 
     vec< vec< pair<int,int> > >& index, const vecbasevector& genome2, const int LG, 
     const VecIntPairVec& Glocs, const vec< digraphVE<int,int> >& H,
     const vec< digraphVE<int,int> >& Hrc, const vec<int>& uperfect,
     const Bool VERBOSE );

void ExtendAlignment( const vec<int>& x1, const vec<int>& x2,
     const vec<int>& nkmers, const int p1, const int p2, const int max_delta,
     vec< pair<int,int> >& a );

void GetAligns( const int id1, const vec<int>& x1, 
     const vec<int>& nkmers, const vec< vec<int> >& ucores,
     const vec< vec< pair<int,int> > >& cindex,
     vec< pair<ualign, int> >& aligns_ids );

void GetAligns( const int id1, const int p1, const vec<int>& x1, 
     const vec<int>& nkmers, const vec< vec<int> >& ucores,
     const vec< vec< pair<int,int> > >& cindex,
     vec< pair<ualign, int> >& aligns_ids );

#endif
