///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MapReadsToHyper.  Map reads to a HyperBasevector.  Crappy interface,
// insufficiently robust for production code, and inadequately documented, but
// potentially useful.

#ifndef MAP_READS_TO_HYPER_H
#define MAP_READS_TO_HYPER_H

#include "Basevector.h"
#include "CoreTools.h"
#include "Qualvector.h"
#include "paths/HyperBasevector.h"

void MapReadsToHyper( const vecbasevector& bases, const vecqualvector& quals, 
     const HyperBasevector& hb, vec< vec< vec<int> > >& rpaths, 
     vec< vec< pair<int,int> > >& aligns, vec< vec<int> >& sum, const int minq );

#endif
