///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file DigraphVTripleIIB.h
 * \author tsharpe
 * \date Jul 6, 2012
 *
 * \brief
 */
#ifndef PATH_DIGRAPHVTRIPLEIIB_H_
#define PATH_DIGRAPHVTRIPLEIIB_H_

#include "graph/Digraph.h"
#include "STLExtensions.h"

typedef triple<int,int,Bool> tripIIB;
typedef digraphV< tripIIB > hyperIIB;
extern template class digraphV< tripIIB >;

#endif /* PATH_DIGRAPHVTRIPLEIIB_H_ */
