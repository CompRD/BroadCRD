/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////



/// Transform one fasta file into vecfasta format, save names separately.
/// 
/// \file Fasta2Vecfasta.cc
///
/// names are saved as a vecString in OUT.names



#include "Basevector.h"
#include "MainTools.h"
#include "FastaFileset.h"
#include "Fastavector.h"

int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String(IN);
  CommandArgument_String_OrDefault(OUT,"");

  CommandArgument_Bool_OrDefault_Doc(NAMES, True,
		     "Whether to save the original sequence names as found in "
				     "fasta file into fastb.names file");
  EndCommandArguments;

  BaseVecVec reads;
  vecfastavector vec_reads;
  vecString readNames;

  if (OUT.empty()) OUT=IN.SafeBefore(".fa") + ".vecfasta";

  FastFetchReads(reads, &readNames, IN);


  for (size_t i = 0; i < reads.size(); i++) 
    vec_reads.push_back_reserve( fastavector( reads[i] ) );
  
  vec_reads.WriteAll(OUT);
  if (NAMES) readNames.WriteAll(OUT + ".names");
  return 0;
}
