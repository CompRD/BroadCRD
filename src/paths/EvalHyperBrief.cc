/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// EvalHyperBrief: print out info about HyperKmerPath component sizes.

#include "MainTools.h"
#include "math/Functions.h"
#include "paths/EvalUtils.h"
#include "paths/HyperKmerPath.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String(SUBDIR);
     CommandArgument_String_OrDefault(WRUN, "run");
     CommandArgument_String_OrDefault(HYPER, "hyper");
     CommandArgument_Int_Doc( K, "kmer size (must match the ones used elsewhere)." );
     EndCommandArguments;

     // Set up directories.

     String data_dir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
     String wdata_dir = sub_dir;
     String wrun_dir = sub_dir + "/" + WRUN;

     // Load HyperKmerPath.

     HyperKmerPath h( sub_dir + "/" + HYPER );
     ForceAssertEq( K, h.K( ) );

     // Compute stats.

     vec<int> component_sizes;
     ComputeComponentSizes( h, component_sizes );
     Sort(component_sizes);
     for ( int i = 0; i < component_sizes.isize( ); i++ )
          component_sizes[i] += K - 1;
     cout << component_sizes.size( ) << " components of N50 size "
          << N50(component_sizes) << " and total size "
          << BigSum(component_sizes) << endl;    }
