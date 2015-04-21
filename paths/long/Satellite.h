///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef SATELLITE_H
#define SATELLITE_H

#include "CoreTools.h"

// targets: list, allowed values are alpha, two, ebv, all

void CleanSatelliteReads( const String& TMP, String targets,
     const double max_alpha_score, const Bool verbose = False );

#endif
