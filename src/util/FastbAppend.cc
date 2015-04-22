// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology
/// Append one fastb file another.
/// 
/// \file FastbAppend.cc
///

#include "Basevector.h"
#include "MainTools.h"

int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String(IN1);
  CommandArgument_String(IN2);
  CommandArgument_String(OUT);
  EndCommandArguments;

  vecbasevector seqs1( IN1 );
  vecbasevector seqs2( IN2 );

  seqs1.Append(seqs2);

  seqs1.WriteAll( OUT );

  longlong totalBases = 0;
  for ( size_t i = 0; i < seqs1.size(); ++i )
    totalBases += seqs1[i].size();

  int totalSeqs = seqs1.size();

  PRINT2( totalSeqs, totalBases );

  return 0;
}
