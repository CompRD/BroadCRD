 ///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/**
   Program: Reads2Unibases

   Computes unibases from a set of reads

   INPUT FILES:
     reads.fastb

   OUTPUT FILES:
     unibases.fastb
     
   @file

*/


#include "MainTools.h"
#include "paths/UnibaseUtils.h"
#include "VecUtilities.h"
#include "paths/Unipath.h"
#include "paths/ReadsToPathsCoreX.h"


int main( int argc, char *argv[] ){
  RunTime( );
  
  BeginCommandArguments;

  CommandArgument_String(READS_IN);
  CommandArgument_String_OrDefault(UNIBASES, "unibases");
  CommandArgument_Int(K);
  CommandArgument_UnsignedInt_OrDefault(N_THREADS, 8);

  EndCommandArguments;
  

  
  
  
  vecbasevector reads(READS_IN);
  
  cout << Date() << ": pathing" << endl;
  vecKmerPath paths, paths_rc;
  vec<big_tagged_rpint> pathsdb;
  ReadsToPathsCoreY( reads, K, paths, paths_rc, pathsdb, "", N_THREADS);
  
  cout << Date() << ": unipathing" << endl;
  vecKmerPath unipaths;
  vec<big_tagged_rpint> unipathsdb;
  Unipath( paths, paths_rc, pathsdb, unipaths, unipathsdb, True);
  
  cout << Date() << ": building unipath adjacency graph" << endl;
  digraph unigraph;
  BuildUnipathAdjacencyGraph( paths, paths_rc, pathsdb, unipaths, unipathsdb, unigraph );
  BinaryWriter::writeFile(  READS_IN.SafeBefore(".fastb") + ".unipath_adjgraph"+ ".k" + ToString(K), unigraph );

  cout << Date() << ": building KmerBaseBroker (needed for unibases)" << endl;
  KmerBaseBrokerBig kbb( K, paths, paths_rc, pathsdb, reads );
      
  cout << Date() << ": building unibases" << endl;
  vecbasevector unibases;
  unibases.reserve(unipaths.size());
  for ( size_t i = 0; i < unipaths.size(); i++ ){
    unibases.push_back( kbb.Seq( unipaths[i] ) );
    basevector unib = unibases.back(), rcunib = unibases.back();
    rcunib.ReverseComplement();
    if ( unib == rcunib ){
      cout << ">" << i << " palindromic. l = " << unib.size() << "\n" << ToString(unib) << endl; 
    }
  }
  
  PRINT( unibases.size() );

  unibases.WriteAll( READS_IN.SafeBefore(".fastb") + "." + UNIBASES + ".k" + ToString(K) + ".fastb");

  cout << Date() << ": building unibase overlap graph at " << K -1 << endl;
  digraph unibgraph;
  BuildUnibaseAdjacencyGraph( unibases, unibgraph, K );
  BinaryWriter::writeFile(  READS_IN.SafeBefore(".fastb") + "." + UNIBASES + ".k" + ToString(K) + ".ovrlpgraph"+ ".k" + ToString(K), unibgraph );

  cout << "Reads2Unibases done!" << endl;
}

