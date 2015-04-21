/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// \file RandomBases.cc.  Generate one or more sequences of N random bases.

// Defaults are set up so that the behavior is unchanged from the
// previous version, which only generated one sequence in fasta format
// to stdout, if none of the new optional arguments are given.

// Parameters:
// - N : length of the sequence(s)
// - NSEQS : number of sequences; default 1.
// - OUT : Name of fastb file to write.  Default "" yields fasta file to stdout.
// - SEED : seed for srandomx.  Default 0 means no call to srandomx.

#include "MainTools.h"
#include "Basevector.h"
#include "random/Random.h"

int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_Int(N);
  CommandArgument_UnsignedInt_OrDefault(NSEQS, 1);
  CommandArgument_String_OrDefault(OUT, "");
  CommandArgument_Int_OrDefault(SEED, 0);
  EndCommandArguments;

  if (SEED!=0) srandomx(SEED);
  
  basevector b(N);
  vecbasevector seqs;
  if (!OUT.empty())
    seqs.Reserve(N*NSEQS, NSEQS);

  for (unsigned int s = 0; s<NSEQS; ++s) {
    for ( int i = 0; i < N; i++ )
      b.Set( i, randomx( ) % 4 );
    if (OUT.empty()) {
      b.Print( cout, ToString(N) + "_random_bases"
	       + ((NSEQS>1) ? ("_sequence_" + ToString(s)) : String("")) );
    } else {
      seqs.push_back(b);
    }
  }

  if (!OUT.empty())
    seqs.WriteAll(OUT);

  return 0;
}
