///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "graph/Digraph.h"
#include "math/NStatsTools.h"
#include "feudal/BinaryStream.h"

/**
 * SimpleDigraphStats
 *
 * Print basic stats for a digraph object.
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( DIGRAPH );
  EndCommandArguments;

  cout << Date( ) << ": loading graph\n" << endl;
  digraph graph; 
  BinaryReader::readFile( DIGRAPH, &graph );

  // Determine graph components

  vec< vec<int> > comps;
  graph.Components( comps );

  size_t disconnected = 0;
  size_t connected = 0;
  vec<int> compVertexCount(comps.size());   // number of vertices in component

  for (size_t i = 0; i < comps.size(); i++) {
    compVertexCount[i] = comps[i].size();
    if (compVertexCount[i] == 1) 
      disconnected++;
    else
      connected += compVertexCount[i];
  }

  size_t totalVertexCount = disconnected + connected;
  cout << "Graph consists of " << comps.size() << " components" << endl;
  cout << "Disconnected vertices : " << disconnected << endl;
  cout << "Connected vertices    : " << connected << endl;
  cout << endl;

  cout << "N-STATS FOR THE COUNT OF VERTICES IN EACH COMPONENT\n" << endl;
  PrintBasicNStats( "# of vertices", compVertexCount, cout );
  cout << endl;
  
  ReverseSort(compVertexCount);
  cout << "10 largest (by vertex count) graph components:" << endl;
  for ( size_t i = 0; i < Min(compVertexCount.size(), (size_t) 10); i++ ) {
    double percSize = 100.0 * (double) compVertexCount[i] / (double) totalVertexCount;
    cout << "\t" << compVertexCount[i] << " vertices,  " << percSize << " % of all vertices" << endl;
  }


  cout << Date( ) << ": done" << endl;
}
