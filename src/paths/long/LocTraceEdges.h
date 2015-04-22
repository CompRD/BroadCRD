///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//


#ifndef LOCTRACEEDGES_H_
#define LOCTRACEEDGES_H_

#include "paths/long/EvalByReads.h"

vec<vec<read_place> > LocTraceEdges( const HyperBasevector& hb_fw, const vecbasevector& bases, const vecqualvector& quals, const double delta_qsum, const int64_t L, bool debug =false);


#endif /* LOCTRACEEDGES_H_ */
