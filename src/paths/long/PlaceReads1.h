///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PLACE_READS1_H
#define PLACE_READS1_H

#include "Basevector.h"
#include "CoreTools.h"
#include "Qualvector.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

void PlaceReads1( const HyperBasevector& hb, const vecbasevector& bases,
     const vecqualvector& quals, ReadPathVec& paths );

#endif
