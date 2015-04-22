// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

#include "MainTools.h"
#include "graphics/PartitionedLineGraph.h"
#include "STLExtensions.h"
#include <algorithm>
#include <numeric>
#include <iostream>


int main( int argc, char** argv )
{
  RunTime();
  BeginCommandArguments;
  CommandArgument_String( title );
  CommandArgument_String( values );
  EndCommandArguments;

  ifstream valueStrm( values.c_str() );

  vec<float> chunks;
  float chunk;
  while ( valueStrm >> chunk )
    chunks.push_back( chunk );

  sort( chunks.rbegin(), chunks.rend() );

  float total = accumulate( chunks.begin(), chunks.end(), 0.0 );

  using namespace ns_psplot;

  int height = 5, width = 500;
  vec<freetext> labels;
  labels.push_back( freetext(string(title), black, 16) );
  labels.push_back( freetext(string("sorted by size from left to right"), black, 10) );

  partitioned_line_graph plg = PartitionedLinePlot( chunks, total, 
						    height, width,
						    labels,
						    true );
  
  ofstream test_out( "test.eps" );
  test_out << plg << endl;
  test_out.close();

  System( "gv test.eps" );

  cout << "Output in test.eps." << endl;

  exit(0);
}
