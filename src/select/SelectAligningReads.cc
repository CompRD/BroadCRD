/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

///Select NAMES of reads from FASTA file that align to the lookup table LOOKUP.
/// \file SelectAligningReads.cc
///
/// Parameters:
/// - FASTA: input fasta file of reads to select from.
/// - LOOKUP: lookup table 
/// - NAMES: output names file with names of selected reads. If empty,
/// it is constructed from the basenames of FASTA and LOOKUP. This file should
/// then be fed into MakeFastaSubset.

#include "MainTools.h"
#include "Basevector.h"
#include "FastaFileset.h"
#include "lookup/LookAlign.h"
#include "lookup/LookAlignSort.h"

int main(int argc, char ** argv) {

  BeginCommandArguments;

  CommandArgument_String(FASTA);
  CommandArgument_String_OrDefault(NAMES,"");
  CommandArgument_String(LOOKUP);

  EndCommandArguments;

  if (NAMES.empty()) {
    NAMES = FASTA.SafeBefore("fa") + "match." 
      + Basename(LOOKUP).SafeBefore(".lookup") + ".names";
  }

  vecbasevector reads;
  vecString names;
  vecqualvector readsq;

  if (FASTA.Contains(".gz")) {
    String command = "zcat " + FASTA;
    FASTA=FASTA.SafeBefore(".gz");
    command = command + " > " + FASTA;
    SystemSucceed(command);
  }

  FastFetchReads(reads, &names, FASTA);

  temp_file qltoutname("SelectAlignRds_qltout_XXXXXX");
  
  cout << "Aligning " << reads.size() << " reads to " << LOOKUP << endl;
  SystemSucceed("QueryLookupTable SEQS=" + FASTA + " L=" + LOOKUP
		+ " K=12 MF=1000:10000:100000 END_STRETCH=2 MO=10 "
		"PARSEABLE=True VISUAL=True VISUAL_ABBR=False KB=10 > " 
		+ qltoutname);

  cout << "Loading read alignment data from " << qltoutname << endl;
  // TODO: potentially dangerous truncation of index by readNumbers
  vec<int> readNumbers(reads.size());
  for (size_t i=0; i != reads.size(); ++i) readNumbers[i]=i;
  vec<look_align_plus> aligns;
  vec<vec<int> > alignIndices;
  LoadLookAlignPlus(qltoutname, aligns);
  GetAlignIndices(readNumbers, aligns, alignIndices);
  SortAlignIndices(aligns, alignIndices);
  cout << "Read " << aligns.size() << " alignments." << endl;

   
  //open the output file.
  Ofstream(namestream, NAMES);

  int found=0;
  vec<String> goodnames;
  for (size_t i=0; i != reads.size(); ++i) {
    if (!alignIndices[i].empty()) {
      goodnames.push_back(names[i]);
      ++found;
    }
  }
  namestream << goodnames;
  
  cout << "Put " << found << " names out of " 
       << reads.size()  << " total reads into " << NAMES << endl;
  return 0;
}
