/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   Program: BuildUnipathAdjGraph

   Builds adjacency graph from unipaths and read paths information or directly from
   unibases (using K-1 mer overlap).

   @file
*/

#include "Alignment.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "MainTools.h"
#include "paths/ImproperMerge.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"
#include "paths/Unipath.h"
#include "paths/simulation/Placement.h"
#include "paths/UnibaseUtils.h"

int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandDoc("Build unipath (or unibase) adjacency graph.");
  
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_Int(K);
  CommandArgument_String_OrDefault(READS, "reads");
  CommandArgument_String_OrDefault(PATHS, "paths");
  CommandArgument_String_OrDefault_Doc(UNIPATHS, "unipaths",
	 "Output unipaths name - READS.UNIPATHS.kK");
  CommandArgument_String_OrDefault_Doc(UNIBASES, "unibases",
	 "Output unibases name - READS.UNIBASES.kK");
  CommandArgument_String_OrDefault_Doc(UNIPATH_GRAPH, "adjgraph",
	 "Output adjacency graph name - READS.UNIPATHS.UNIPATH_GRAPH.kK");
  CommandArgument_String_OrDefault_Doc(UNIBASE_GRAPH, "ovrlp_adjgraph",
	 "Output adjacency graph name - READS.UNIBASES.UNIBASE_GRAPH.kK");
  CommandArgument_Bool_OrDefault_Doc(BUILD_UNIPATH_GRAPH, True,
	 "Build and write the unipath adjacency graph. Otherwise build and write unibase adjacency graph");

  EndCommandArguments;
  
  String data_dir = PRE + "/" + DATA;
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  
  String KS = ToString(K);
  

  digraph unigraph;

  // Read in data.  
 
  if ( BUILD_UNIPATH_GRAPH ){
 
    cout << Date( ) << ": loading pathing data" << endl;
    vecKmerPath paths( run_dir + "/" + READS + "." + PATHS + ".k" + ToString(K) );
    vecKmerPath paths_rc( run_dir + "/" + READS + "." + PATHS + "_rc.k" + ToString(K) );
    vecKmerPath unipaths( run_dir + "/" + READS + "." + UNIPATHS + ".k" + ToString(K) );
    BREAD2( run_dir + "/" + READS + "." + UNIPATHS + "db.k" + ToString(K), 
	    vec<tagged_rpint>, unipathsdb );
    BREAD2( run_dir + "/" + READS + "." + PATHS + "db.k" + KS, vec<tagged_rpint>, pathsdb );

    // Build unipath adjacency graph.
    cout << Date() << ": building unipath adjacency graph" << endl;
    BuildUnipathAdjacencyGraph( paths, paths_rc, pathsdb, unipaths, unipathsdb, unigraph );
  }else{
    
    cout << Date( ) << ": loading unibases" << endl;
    vecbasevector unibases( run_dir + "/" + READS + "." + UNIBASES + ".k" + ToString(K) );
    int nuni = unibases.size();
    cout << Date() << ": building unibase adjacency graph" << endl;
    BuildUnibaseAdjacencyGraph( unibases, unigraph, K );
  }
  BinaryWriter::writeFile( run_dir + "/" + READS + "." + UNIBASES + "." + UNIBASE_GRAPH + ".k" + KS, unigraph );

  cout << Date( ) << ": Done with BuildUnipathAdjGraph!" << endl;
  return 0;
}
