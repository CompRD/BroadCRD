/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef MATCHES_PERF_H
#define MATCHES_PERF_H

#include "Vec.h"

/**
 * MatchesPerf
 *
 * Return true iff bases1 and bases2 match perfectly. Warning: if
 * bases2 is shorter than bases1, then it will test for perfect match
 * by "attaching" bases2 to bases1 in the two possible ways (if
 * instead bases1 is shorter than bases2, then MatchesPerf return
 * false).
 */
bool MatchesPerf( const vec<char> &bases1,
		  const vec<char> &bases2,
		  ostream *plog = 0 );

#endif
