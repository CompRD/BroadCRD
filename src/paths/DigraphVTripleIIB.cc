///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file DigraphVTripleIIB.cc
 * \author tsharpe
 * \date Jul 6, 2012
 *
 * \brief
 */

#include "paths/DigraphVTripleIIB.h"
#include "graph/DigraphTemplate.h"
template digraphV<tripIIB>::digraphV(const vec<vec<int> >&, const vec<vec<int> >&, const vec<tripIIB>&);
template int digraphV<tripIIB>::AddVertex(const tripIIB&);
template void digraphV<tripIIB>::DeleteVertex(const int);
template void digraphV<tripIIB>::DeleteVertices(const vec<int>&);
template void digraphV<tripIIB>::Initialize(vec<vec<int> > const&, vec<vec<int> > const&, vec<triple<int, int, unsigned char> > const&);
template const tripIIB& digraphV<tripIIB>::Vert(int) const;
