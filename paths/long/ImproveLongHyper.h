///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef IMPROVE_LONG_HYPER_H
#define IMPROVE_LONG_HYPER_H

#include "Basevector.h"
#include "CoreTools.h"
#include "Qualvector.h"
#include "efasta/EfastaTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/DataSpec.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/SupportedHyperBasevector.h"

void ImproveLongHyper( const String& SAMPLE, const String& X, 
     SupportedHyperBasevector& shb, const String& TMP, const String& READS,
     const long_data_spec& spec, const long_heuristics& heur, unsigned nThreads, 
     const long_logging_control& log_control, const long_logging& logc, 
     const int start_step, bool useOldLRPMethod );

#endif
