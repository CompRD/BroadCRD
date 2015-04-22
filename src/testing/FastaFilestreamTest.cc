// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
//


#include "FastaFilestream.h"
#include "system/RunTime.h"
#include <iostream>

int main( int argc, char** argv )
{

  RunTime();

  if ( argc != 2 )
  {
    cerr << "Specify a filename, please." << endl;
    exit( -1 ) ;
  }

  String filename( argv[1] );

  FastaQualityFilestream fasta_test( filename, new FullNameParser );

  vec<int> read_indices;
  read_indices.push_back( 1 );
  read_indices.push_back( 4 );
  sort( read_indices.begin(), read_indices.end() );

  vecString names;
  vecqualvector quals;

  fasta_test.parseSubset( read_indices, names, quals );

  for ( vecqvec::size_type i = 0; i < quals.size(); ++i )
  {
    cout << names[i] << endl;
    for ( qvec::size_type j = 0; j < quals[i].size(); ++j )
      cout << (int) quals[i][j];
    cout << endl;
  }

  exit(0);
}
