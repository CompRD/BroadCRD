// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research


#include "MainTools.h"
#include "STLExtensions.h"
#include <algorithm>
#include <numeric>

int main( int argc, char** argv )
{
  RunTime();

  vector<int> V(10);
  iota( V.begin(), V.end(), 1 );
  
  copy( V.begin(), V.end(), ostream_iterator<int>( cout, " " ) );
  cout << endl;

  vector<int> Vmin;
  copy( V.begin(), V.end(), back_inserter(Vmin) );

  transform( Vmin.begin(), Vmin.end(), Vmin.begin(), bind2nd(minimum<int>(), 5) );

  copy( Vmin.begin(), Vmin.end(), ostream_iterator<int>( cout, " " ) );
  cout << endl;

  vector<int> Vmax;
  copy( V.begin(), V.end(), back_inserter(Vmax) );

  transform( Vmax.begin(), Vmax.end(), Vmax.begin(), bind2nd(maximum<int>(), 5) );

  copy( Vmax.begin(), Vmax.end(), ostream_iterator<int>( cout, " " ) );
  cout << endl;

  exit(0);
}
