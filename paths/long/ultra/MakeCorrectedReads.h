///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef MAKE_CORRECTED_READS_H
#define MAKE_CORRECTED_READS_H

#include "Basevector.h"
#include "CoreTools.h"
#include "efasta/EfastaTools.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/Simulation.h"

void MakeCorrectedReads( const long_sim& sim, const vecbasevector& reads,
     const String& RID, const unsigned int NUM_THREADS, const long_heuristics& heur,
     const long_logging_control& log_control, const long_logging& logc, 
     const ref_data& ref, VecEFasta& corrected, vec<int>& cid, String const& TMP );

#endif
