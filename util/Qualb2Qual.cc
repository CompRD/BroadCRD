/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/// Transform qualb file IN into fasta format.
///
/// \file Qualb2Qual.cc
///
/// Default for OUT is IN's prefix + .qual
/// If GENERATE_NAMES is set names must be in file IN.names or IN.names.gz.

#include "Qualvector.h"
#include "MainTools.h"
#include "FetchReads.h"

int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String(IN);
  CommandArgument_String_OrDefault(OUT,"");
  CommandArgument_Bool_OrDefault(GENERATE_NAMES,True);
  EndCommandArguments;

  if (OUT.empty()) OUT= IN.SafeBefore(".qu")+".qual";

  vecqualvector reads;
  reads.ReadAll(IN);
  vecqvec::size_type N=reads.size();

  vecString readNames;
  if (GENERATE_NAMES) {
    for (vecqvec::size_type i=0; i<N; i++)
      readNames.push_back_reserve(ToString(i));
  } else {
    readNames.ReadAll(IN+".names");
    ForceAssertEq(N, readNames.size());
  }

  Ofstream(os,OUT);
  for (vecqvec::size_type i = 0; i != N; ++i) {
    Print(os, reads[i], readNames[i]);
  }

  return 0;
}
