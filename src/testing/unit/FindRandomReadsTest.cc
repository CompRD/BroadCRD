/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


#include  "random/FindRandomReads.h"
#include "MainTools.h"

int main(int argc, char ** argv) {
  srand48(time(NULL));
  vec<int> sizes;
  sizes.push_back(10000,1000,100,10);
  vec<pair<int, int> > positions;
  FindRandomReads::Positions(10000,sizes,positions);
  vec<int> counts(sizes.size(),0);
  for (int i=0; i != positions.isize(); ++i) {
    counts[positions[i].first]++;
  }
  PRINT(counts);
  return 0;
}
