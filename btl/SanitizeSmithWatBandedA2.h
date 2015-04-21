///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef BTL__SANITIZE_SMIT_WAT_BANDED_A2__H
#define BTL__SANITIZE_SMIT_WAT_BANDED_A2__H

#include "Alignment.h"

/**
 * SanitizeSmithWatBandedA2
 *
 * Walk around a known bug in SmithWatBandedA2's output, which causes
 * the last block to have size 0 in some cases. This in turn would
 * cause some other code (eg LoadLookAligns) to crash.
 *
 * Note that the align will be replaced. The code returns false if any
 * other (unfixable) inconsistency is found.
 */
bool SanitizeSmithWatBandedA2( align &al );

#endif
