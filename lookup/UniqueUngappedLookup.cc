/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


/** Find good unique alignments with no indels quickly.  Intended to
    replace ImperfectLookupTable.  

\file UniqueUngappedLookup.cc

Use LookupQuery() and UniqueGlobalUngappedAlignsFromHits to find alignments. Save two files: 

- OUT_PREFIX.qltout contains the unique alignments in parseable 
and readable brief format.
- OUT_PREFIX.minErrors.txt: one line per read, with the number of errors 
in the best alignment for each read, or -1 if there are no 
alignments.

Parameters:
- SEQS: fasta or fastb file
- LOOKUP_TABLE, L: file name for lookup table. Note the table's K size is used for the lookup.
- OUT_PREFIX, O: prefix for output files, by default same as fasta file.
- OUT_SUFFIX: (default uul.qltout) suffix for output alignment file.
- FWRC (default True): align in both direction, if False FW only.
- ERR_DIFF (default 2): if the second best alignment has ERR_DIFF or fewer
errors more than the best alignment, do not save any alignments for this
read (but do save the number of errors in the best alignment).
- MAX_ERRS (default 4): if the best alignment has more than MAX_ERRS
errors, do not save any alignments for this read (but do save the
number of errors in the best alignment).
- CHUNK (default 500000): use this many reads at a time. This allows us
to limit the memory demands, which is critical when running the 
pipeline on the blades.

*/

#include "MainTools.h"
#include "lookup/HitReceiver.h"
#include "FastaFileset.h"

int main( int argc, char *argv[] )
{  
  RunTime();

  BeginCommandArguments;

  CommandArgument_String(SEQS);
  CommandArgument_String_Abbr(LOOKUP_TABLE, L);
  CommandArgument_String_Abbr(OUT_PREFIX, O);
  CommandArgument_String_OrDefault(OUT_SUFFIX, ".uul.qltout");
  CommandArgument_Bool_OrDefault(FWRC, True);
  CommandArgument_Bool_OrDefault(WRITE_MINALIGNS, True);
  CommandArgument_Int_OrDefault(ERR_DIFF,2);
  CommandArgument_Int_OrDefault(MAX_ERRS,4);
  CommandArgument_UnsignedInt_OrDefault(MAX_FREQ,0);
  CommandArgument_UnsignedInt_OrDefault(CHUNK,500000);
  EndCommandArguments;

  if (OUT_PREFIX.empty()) {
    OUT_PREFIX = SEQS.SafeBefore(".fa");
  }

  vecbasevector seqs;

  if (SEQS.Contains(".fastb")) {
    seqs.ReadAll(SEQS);
  } else {
    FastFetchReads(seqs, 0, SEQS);
  }
  
  Ofstream(out, OUT_PREFIX + OUT_SUFFIX);
  PrintCommandPretty(out);

  ofstream * errStream = 0;
  if (WRITE_MINALIGNS) {     
    errStream = new ofstream( (OUT_PREFIX + ".minAlignErrors.txt").c_str());
  }

  cout << Date() << ": input read.\n";
  for (size_t c=0; c != (seqs.size() / CHUNK) + 1; ++c) {
    // TODO: potentially dangerous truncation of index by firstRead/lastRead
    size_t firstRead = c*CHUNK;
    size_t lastRead = min((c+1)*CHUNK, seqs.size());
    PRINT2(firstRead, lastRead);
    double startTime = WallClockTime();
    vec<look_align> aligns; ///< Results from lookup table, indexed by query seq
    vec<int> min_errors;   ///< Results from lookup table, indexed by query seq
    vec<int> next_best_errors;   ///< Results from lookup table, indexed by query seq
    { // Read in lookup table and get hits against it
      lookup_table look(LOOKUP_TABLE);
      UniqueGlobalUngappedHitReceiver align(look, seqs, aligns, min_errors, next_best_errors,
					    ERR_DIFF, MAX_ERRS);
      look.FindHits(seqs, align, MAX_FREQ, (FWRC ? 2 : 1),
		    firstRead, lastRead);
      cout << "Hits took " << TimeSince(startTime) << endl;
    }

    // Now done with lookup table; write out results
    if (WRITE_MINALIGNS) {
      for (int i=0; i != min_errors.isize(); ++i) {
	if ((unsigned int)min_errors[i] == infinitely_many) min_errors[i] = -1;
      }
      copy(min_errors.begin(), min_errors.end(), 
	   ostream_iterator<int>(*errStream,"\n"));
    }

    for (int i=0; i != aligns.isize(); ++i) {
      if ( 0 <= min_errors[i] && min_errors[i] <= MAX_ERRS &&
	   (next_best_errors[i] - min_errors[i]) > ERR_DIFF ) {
	aligns[i].PrintParseable(out);
	aligns[i].PrintReadableBrief(out);
      }
    }
    out.flush();
  }
  delete errStream;

  return 0;
}
