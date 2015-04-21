/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// Program: FindKmerFrequencies
//
// See Also: FindStrongKmers, FindGenomicKmers
//
// Find all the k-mers in a given fastb file and its reverse complement.
// Output these kmers and their multiplicities, capped at 65535, optionally 
// excluding unique kmers.


#include "MainTools.h"
#include "kmer_freq/WriteKmerFrequencies.h"
#include "kmers/KmerShape.h"

int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_KShape(K);
  CommandArgument_String_Doc(SEQS_IN,
    "Input fastb file from which to generate the kmer frequency table");
  CommandArgument_String_OrDefault_Doc(KMER_FILE_OUT, "",
    "Ouput kmer freq table filename - otherwise append .freq_table.kn");
  CommandArgument_Bool_OrDefault_Doc(INCLUDE_UNIQUE, True,
    "Include uniquely occuring kmers in frequency table");
  EndCommandArguments;

  String kmerFile;
  
  if (KMER_FILE_OUT == "")
    kmerFile = SEQS_IN.SafeBefore(".fastb") + ".freq_table.k" + ToString(K);
  
  vecbasevector bases(SEQS_IN);

  cout << kmerFile << endl;
#define CASE(_KSHAPE) \
  WriteKmerFrequencies<_KSHAPE>(bases, kmerFile , INCLUDE_UNIQUE )
  
  DISPATCH_ON_KSHAPE(K, CASE);

}

