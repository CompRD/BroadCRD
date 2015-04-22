/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "FastaFileset.h"
#include "VecUtilities.h"
#include <cstdlib>

/** Pick LINES lines randomly from IN without replacement NUM_FILES times.
And put these sets of lines into NUM_FILES separate files.

\file SubsetFastbFile.cc

*/

int main( int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String(IN); 
  CommandArgument_Int_OrDefault(SIZE,100000); 
  CommandArgument_Int_OrDefault(SEED,177777); 
  CommandArgument_String(OUT); 
  EndCommandArguments;

  vecbasevector reads(IN);

    vec<size_t> randoms(reads.size());
    vec<size_t> perm(reads.size());
    for ( size_t i = 0; i < reads.size(); i++)
    {
        randoms[i] = rand();
        perm[i] = i;
    }

    SortSync(randoms, perm);

    vecbasevector subset(SIZE);
    for (int i = 0; i < SIZE; i++)
    {
        subset[i] = reads[perm[i]];
    }
    subset.WriteAll(OUT);

  return 0;
}

