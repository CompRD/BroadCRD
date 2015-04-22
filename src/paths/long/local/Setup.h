///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LOCAL_SETUP_H
#define LOCAL_SETUP_H

// MakeDepend: private

#include "CoreTools.h"
#include "paths/long/DataSpec.h"
#include "paths/long/Logging.h"

// Specification of a Picard dataset.  If special is specified, it is taken
// to be the path of a bam file, and the other arguments are ignored.

class picard_spec {

     public:

     picard_spec( ) { }
     picard_spec( const Bool direct, const String& flowcell, const String& picard_run,
          const int lane, const String& lib, const String& special = "" )
          : direct(direct), flowcell(flowcell), picard_run(picard_run), lane(lane), 
          lib(lib), special(special) { }

     Bool direct;
     String flowcell;
     String picard_run;
     int lane;
     String lib;
     String special;
};

void AddTo( const Bool direct, vec<picard_spec>& specs, const String& flowcell, 
     const String& picard_run, const int lane, const String& lib );

void AddTo( const Bool direct, vec<picard_spec>& specs, const String& flowcell, 
     const String& picard_run, const vec<int>& lanes, const String& lib );

void AddTo( const Bool direct, vec<picard_spec>& specs, const String& flowcell, 
     const String& picard_run, const int lane, const vec<String>& libs );

inline int ExtendRegionBy( ) { return 400; }

void SetupIlluminaData( 
     // inputs:
     const long_data_spec& spec, const String& SAMPLE, const String& DATASET, 
     const vec<String>& regions, const double COVERAGE_CAP, const Bool HPOOL_ORIG,
     // outputs:
     vec<picard_spec>& specs, vec<String>& chrlist, String& region,
     // input and output:
     double& SELECT_FRAC,
     // logging:
     const long_logging& logc );

void Dexterize( const String& var, String& value, String& SAMPLE,
     String& X, double& GENOME_SUB_PERCENT, String& TMP, String& READS );

void DexterizeAll( String& IN_SHBV, String& IN_SHBV_FINAL, String& IN_EFASTA_READS,
     const int START_STEP, String& SAMPLE, String& X, double& GENOME_SUB_PERCENT, 
     String& TMP, String& READS );

#endif
