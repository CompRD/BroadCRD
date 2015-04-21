/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef PADS_RATIO_H
#define PADS_RATIO_H

#include "Vec.h"

/**
 * PadsRatio
 *
 * Return the ratio between pads and total length for a given read. Total
 * length is defined as the unpadded length of "consensus". Convention:
 * if unpadded consensus length is 0, then PadsRatio returns 1.0.
 *
 * Remark: consensus and read must have the same size!
 */
float PadsRatio( const vec<char> &consensus, const vec<char> &read );

#endif
