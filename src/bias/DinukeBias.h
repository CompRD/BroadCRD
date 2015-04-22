/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// DinukeBias.  Given a set of reads, aligned to a reference, print out statistics 
// regarding the "dinucleotide bias" of their start points.

#ifndef DINUKE_BIAS_H
#define DINUKE_BIAS_H

#include "Basevector.h"
#include "CoreTools.h"
#include "lookup/LookAlign.h"

void DinukeBias( const vecbasevector& reads, const vecbasevector& ref,
     const vecbasevector& rref, const vec<look_align>& aligns );

#endif
