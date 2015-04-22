// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology
/// Append one fastb file another.
///
/// \file QualbAppend.cc
///

char const* DOC =
"Appends one qualb file to another and writes the result to disk.  "
"The output file can be the same as either input file.";

#include "Qualvector.h"
#include "MainTools.h"

int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String_Doc(IN1, "The contents of this file will be first.");
  CommandArgument_String_Doc(IN2, "The contents of this file will be second.");
  CommandArgument_String_Doc(OUT, "The result will be written to this file.");
  EndCommandArguments;

  vecqualvector seqs1( IN1 );
  vecqualvector seqs2( IN2 );

  seqs1.Append(seqs2);

  seqs1.WriteAll( OUT );

  longlong totalBases = 0;
  for ( vecqvec::size_type i = 0; i < seqs1.size(); ++i )
    totalBases += seqs1[i].size();

  vecqvec::size_type totalSeqs = seqs1.size();

  PRINT2( totalSeqs, totalBases );

  return 0;
}
