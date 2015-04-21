///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/long/magic/BasicScore.h"

int BasicScore( const basevector& b1, const basevector& b2 )
{    alignment a;
     int best_loc;
     int errs;
     if ( b1.size( ) <= b2.size( ) )
          errs = SmithWatFree( b1, b2, best_loc, a, True, True );
     else errs = SmithWatFree( b2, b1, best_loc, a, True, True );
     return errs;    }
