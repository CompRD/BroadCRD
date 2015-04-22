///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef MAKE_SIM_READS_H
#define MAKE_SIM_READS_H

#include "Basevector.h"
#include "CoreTools.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/Simulation.h"

void MakeSimReads( long_sim sim, const ref_data& ref,
     vecbasevector& reads, vec<ref_loc>& readlocs );

void FilterSimReads( const long_sim& sim,
     const vecbasevector& G, const vec<ref_loc>& readlocs, String& RID );

#endif
