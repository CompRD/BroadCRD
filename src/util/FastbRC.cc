/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/// Reverse complement sequences from a fastb file and save them into a new fastb file.
/// 
/// \file FastbRC.cc
///
/// Default for OUT is IN's prefix + _rc.fastb

#include "Basevector.h"
#include "MainTools.h"
#include "FetchReads.h"

int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String(IN);
  CommandArgument_String_OrDefault(OUT, "");
  EndCommandArguments;

  if (OUT.empty()) OUT=IN.SafeBefore(".fastb")+"_rc.fastb";

  vecbasevector reads;
  reads.ReadAll(IN);

  ReverseComplement(reads);
  reads.WriteAll(OUT);
  
  return 0;
}
