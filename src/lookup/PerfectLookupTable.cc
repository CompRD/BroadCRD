/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/** Wrapper around PerfectLookup(): find perfect alignments with no indels.

\file PerfectLookupTable.cc

Use PerfectLookup() to find alignments. Save two files: 

- OUT_PREFIX.qltout contains the unique alignments in parseable 
and readable brief format.


Parameters:
- SEQS: fasta or fastb file
- K: kmer size, must match lookup table
- LOOKUP_TABLE, L: file name for lookup table
- OUT_PREFIX, O: prefix for output files, by default same as fasta file.
- OUT_SUFFIX: (default qltout) suffix for output alignment file.
- FWRC (default True): align in both direction, if False FW only.


*/

#include "MainTools.h"
#include "lookup/PerfectLookup.h"
#include "FastaFileset.h"

int main( int argc, char *argv[] )
{    
  RunTime();

  BeginCommandArguments;

  CommandArgument_String(SEQS);
  CommandArgument_UnsignedInt(K);
  CommandArgument_String_Abbr(LOOKUP_TABLE, L);
  CommandArgument_String_Abbr_OrDefault(OUT_PREFIX, O, "");
  CommandArgument_String_OrDefault(OUT_SUFFIX, "qltout");
  CommandArgument_Bool_OrDefault(FWRC, True);

  EndCommandArguments;

  if (OUT_PREFIX.empty()) {
    OUT_PREFIX = SEQS.SafeBefore(".fa");
  }

  vecbasevector seqs;

  if (SEQS.Contains("fasta") ) {
    FastFetchReads(seqs, 0, SEQS);
  } else {
    seqs.ReadAll(SEQS);
  }
  
  if (OUT_PREFIX.empty()) {
    OUT_PREFIX = SEQS.SafeBefore(".fa");
  }
  
  vec<look_align> aligns;
  vec<int> min_errors;
  PerfectLookup(K, seqs, LOOKUP_TABLE, 
		  aligns, 
		  FWRC ? FW_OR_RC : FW);

  Ofstream(pltout, OUT_PREFIX + "." + OUT_SUFFIX);
  for (int i=0; i != aligns.isize(); ++i) {
    aligns[i].PrintParseable(pltout);
    aligns[i].PrintReadableBrief(pltout);
  }

  return 0;
}
