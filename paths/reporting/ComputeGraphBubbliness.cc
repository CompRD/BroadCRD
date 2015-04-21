/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   Program: GomputeGraphBubbliness
   
   @file
*/

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "graph/Digraph.h"
#include "paths/Unipath.h"
#include "feudal/BinaryStream.h"


int main( int argc, char *argv[] ){
  
  RunTime( );
  
  BeginCommandArguments;
  CommandDoc("Compute graph \"bubbliness\".");

  CommandArgument_String_OrDefault_Doc(GRAPH_IN,"",
         "Head of input files.");

  EndCommandArguments;
     


  String graphInFile    = GRAPH_IN;

  
  cout << Date() << ": reading unipath adjacency graph= " + graphInFile << endl;
  digraph AG;
  BinaryReader::readFile( graphInFile, &AG );
  
  cout << Date() << ": building digraphE" << endl;
  digraphE<int> DG( AG );

  cout << Date() << ": computing bubbliness" << endl;

  size_t nvertices = DG.N();
  
  size_t npossible = 0;
  size_t nbubbles  = 0;
  vec<int> targets, ustargets;
  for ( size_t v = 0; v < nvertices; v++ ){
    if ( DG.From(v).size() > 1 ){
      npossible++;
      targets = DG.From(v);
      ustargets = targets;
      UniqueSort( ustargets );
      if ( ustargets.size() < targets.size() ){
	nbubbles++;
      }
    }
  }
  
  double bubbliness  = (double) nbubbles / (double) npossible * 100;
  double branchiness = (double) npossible / (double) nvertices * 100;

  cout <<
    "nvertices - number of vertices in the derived edge graph\n"
    "npossible - number of vertices with at least 2 outgoing edges\n"
    "nbubbles  - number of vertices with 2 outgoing edges converging on the same vertex\n"
    "bubbliness = nbubbles / npossible * 100\n"
    "branchiness = npossible / nvertices * 100\n"
       << endl;

  PRINT5( nvertices, npossible, nbubbles, bubbliness, branchiness );
  
  cout << Date() << ": ComputeGraphBubbliness done!" << endl;
}
