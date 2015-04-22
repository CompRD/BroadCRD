///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef SNORGLE_TOOLS_H
#define SNORGLE_TOOLS_H

#include "Basevector.h"
#include "CoreTools.h"
#include "paths/long/ReadStack.h"
#include "paths/long/large/explore/SnorgleTools.h"

void FindPaths( const readstack& s, vec<basevector>& p, const Bool FP_LOGGING );

#endif
