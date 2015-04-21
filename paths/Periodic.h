///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// IsPeriodic: determine if a basevector is periodic with periodicity between
// 1 and top_period.  For example,
// AAAAAA is periodic with period 1
// ATATATA is periodic with period 2.

#ifndef PERIODIC_H
#define PERIODIC_H

#include "Basevector.h"
#include "CoreTools.h"

Bool IsPeriodic( const basevector& x, const int top_period );

#endif
