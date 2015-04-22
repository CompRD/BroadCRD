/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef APPLY_CORRECTIONS_H
#define APPLY_CORRECTIONS_H

#include "Basevector.h"
#include "Corrector.h"
#include "Qualvector.h"
#include "ReadLocation.h"

/**
 * ApplyCorrections
 *
 * Apply given corrections to a set of contigs. Notice that locs, bases,
 * and quals are both input and output. Quals are capped!
 *
 * snp_info_out: if passed save snp info (human readable)
 */
void ApplyCorrections( const vec<corrector> &fixes,
		       const vec<int> & first_locs,
		       vec<read_location> &locs,
		       vecbasevector &bases,
		       vecqualvector &quals,
		       ostream *snp_info_out = 0 );

#endif
