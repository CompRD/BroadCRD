/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/// /file UnipathGcBias
///
/// Estimate and store the GC bias for the given unipaths.

#include "MainTools.h"

#include "paths/KmerPath.h"
#include "paths/KmerBaseBroker.h"

int main( int argc, char** argv ) {
  
  RunTime();

  BeginCommandArguments;

  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);

  CommandArgument_Int(K);

  CommandArgument_String_OrDefault(BIAS_CURVE,"reads.gc_bias");

  EndCommandArguments;

  String runDir = PRE + "/" + DATA + "/" + RUN;

  vecKmerPath unipaths( runDir + "/reads.unipaths.k" + ToString(K) );
  
  KmerBaseBroker kbb( runDir, K );
  
  vec<double> biasCurve;
  READX( runDir + "/" + BIAS_CURVE + ".k" + ToString(K), biasCurve );
  
  vec<double> unipathBias( unipaths.size() );

  // For each unipath, calculate the average GC bias across all the
  // kmers in the unipath.
  for ( size_t uid = 0; uid < unipaths.size(); uid++ ) {
    basevector bases = kbb.Seq( unipaths[uid] );
    int numKmers = bases.size() - K + 1;
    ForceAssertGe( numKmers, 0 );
    int gc = bases.GcBases( 0, K );
    double biasAve = 0;
    for ( int kstart = 0; kstart < numKmers; ++kstart ) {
      biasAve += biasCurve[gc];
      gc -= IsGC( bases[kstart] );
      if ( kstart + K < (int) bases.size() )
        gc += IsGC( bases[kstart+K] );
    }
    biasAve /= (double) numKmers;
    unipathBias[uid] = biasAve;
    PRINT2( uid, unipathBias[uid] );
  }

  String out = runDir + "/reads.unipaths.gc_bias.k" + ToString(K);
  ofstream outstrm( out.c_str() );
  outstrm << unipathBias;
  outstrm.close();

  return 0;
}
      
