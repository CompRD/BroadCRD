///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "Generates unibases directly from a fastb, without generating paths or "
  "unipaths. Uses minimal memory and is faster than pathing first.";

#include "MainTools.h"
#include "paths/Unipath.h"
#include "paths/ReadsToPathsCoreX.h"
#include "kmers/ReadPather.h"
#include "kmers/ReadPatherDefs.h"

template <unsigned K> 
vecbasevector buildReadPathUnibases(const string& fastb_in,
				    const unsigned NUM_THREADS,
				    const unsigned coverage) {

  cout << Date() << " Building K=" << K << " unipaths " << endl;

  String graphInfoFilename = UnipathGraph<K>::getInfoFilename(fastb_in);
  String countsFile = UnipathGraph<K>::getCountsFilename(graphInfoFilename);
  size_t nKmers = 4 * MasterVec<bvec>::MastervecFileRawCount(fastb_in) / coverage;
  KmerDict<K>* pDict;
  pDict = new KmerDict<K>(nKmers);
  pDict->process(fastb_in,NUM_THREADS);
  UnipathGraph<K>::create(graphInfoFilename,*pDict);
  
  // Generate K+1 unibases

  cout << Date() << " Building K=" << K << " unibases " << endl;

  UnipathGraph<K> graph(graphInfoFilename);
  UnipathEdgeVec const& edges = graph.getAllEdges();

  vecbasevector unibases;
  unibases.reserve(edges.size());
   
  typedef HugeBVec::const_iterator Itr;
  for ( size_t idx = 0; idx < edges.size(); ++idx ) {
    basevector bv;
    UnipathEdge const& edge = edges[idx];
    Itr bItr = graph.getBases(edge.getInitialKmerID());
    Itr bEnd = graph.getBases(edge.getFinalKmerID()) + K;
    bv.append(bItr,bEnd);
    unibases.push_back(bv);
  }

  return unibases;
} 


int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String_Doc(FASTB_IN,
    "Fastb file to generate unibases from");
  CommandArgument_String_OrDefault_Doc(UNIBASES_OUT, "",
    "Unibases filename");
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  CommandArgument_Bool_OrDefault(VALIDATE,False);
  EndCommandArguments;

  // Thread control.
  
  NUM_THREADS = configNumThreads(NUM_THREADS);

  // Hard coded K values for now.

  unsigned const K = 40;       // Output unibase K
  unsigned const K_1 = K + 1;  // Internal K for pathing

  if (UNIBASES_OUT == "")
    UNIBASES_OUT = FASTB_IN.ReplaceExtension(".fastb",".unibases.k" + ToString(K));

  // Build internal unipath graph using ReadPather at K+1

  cout << Date() << " Building K=" << K_1 << " unipaths " << endl;

  VirtualMasterVec<bvec> reads(FASTB_IN);
  size_t nKmers = 4 * reads.sizeSum() / 50; // 50 is a coverage estimate
  UnipathGraph<K_1> graph(reads,VALIDATE,NUM_THREADS,nKmers);
  
  // Generate K+1 unibases

  cout << Date() << " Building K=" << K_1 << " unibases " << endl;

  UnipathEdgeVec const& edges = graph.getAllEdges();

  vecbasevector K_1_unibases;
  K_1_unibases.reserve(edges.size());
   
  typedef HugeBVec::const_iterator Itr;
  for ( size_t idx = 0; idx < edges.size(); ++idx ) {
    EdgeID edgeID(idx);
    bvec bv(graph.getBases(edgeID),graph.getBasesEnd(edgeID));
    K_1_unibases.push_back(bv);
  }

  // Generate K paths and unipaths from K+1 unibases

  cout << Date() << " Pathing K=" << K_1 << " unibases at K=" << K << endl;

  vecKmerPath new_paths, new_pathsrc;
  vec<tagged_rpint> new_pathsdb, new_unipathsdb;
  ReadsToPathsCoreY( K_1_unibases, K, new_paths, new_pathsrc, new_pathsdb,
		     FASTB_IN + ".SeqsToUnipaths", NUM_THREADS);

  cout << Date() << " Building K=" << K << " unipaths" << endl;

  vecKmerPath new_unipaths;
  Unipath( new_paths, new_pathsrc, new_pathsdb, new_unipaths, new_unipathsdb );
  KmerBaseBroker new_kbb( K, new_paths, new_pathsrc, new_pathsdb, K_1_unibases );

  // Generate and write K unibases

  cout << Date() << " Saving K=" << K << " unibases" << endl;

  IncrementalWriter<basevector> writer(UNIBASES_OUT, new_unipaths.size());
  for ( size_t i = 0; i < new_unipaths.size( ); i++ )
    writer.add(new_kbb.Seq( new_unipaths[i] ));
  writer.close();

  
  cout << Date( ) << " Done." << endl;
}


