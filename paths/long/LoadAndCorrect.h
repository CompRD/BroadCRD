///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LOAD_AND_CORRECT_H
#define LOAD_AND_CORRECT_H

#include "Basevector.h"
#include "CoreTools.h"
#include "efasta/EfastaTools.h"
#include "paths/long/DataSpec.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/PairInfo.h"
#include "paths/long/local/Setup.h"

void LoadPacBio( const String& X, const String& TMP, vecbasevector& pb,
     const long_logging_control& log_control, const long_logging& logc,
     const Bool PACBIO_MORE );

const int Fosmid_trim_back = 600;

void LoadAndCorrectIllumina( const String& SAMPLE, const String& READS, 
     const String& DATASET, const String& X, const String& TMP, 
     const long_data_spec& spec, VecEFasta& corrected, vec<int>& cid, 
     vec<pairing_info>& cpartner, const Bool HAVE_PICARD, 
     const long_heuristics& heur, const long_logging_control& log_control, 
     const long_logging& logc, const int NUM_THREADS, const String& EXIT, 
     const double clock, bool useOldLRPMethod, const Bool KEEP_NAMES );

void LoadIlluminaReads( const String& SAMPLE, const vec<picard_spec>& specs,
     const String& region, const String& TMP, const long_logging& logc,
     const Bool USE_PF_ONLY, const Bool KEEP_NAMES );

#endif
