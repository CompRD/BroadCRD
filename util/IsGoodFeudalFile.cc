// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#include "MainTools.h"

int main( int argc, char** argv )
{
  if ( argc != 2 )
  {
    cout << "Expected one argument." << endl;
    exit( 1 );
  }

  if ( IsGoodFeudalFile( argv[1] ) )
    cout << "yes" << endl;
  else
    cout << "no" << endl;
  
  exit( 0 );
}
