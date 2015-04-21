/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/** Print out how many bases have coverage 0,1,... n
 \file CoverageHistFromAligns.cc

This gives a quick readout of coverage for human consumption and
estimation of randomness. Uses CoverageAnalyzer.

   - QLTOUT: list of text look_align files, separated by DELIMs.
   - REF: fastb file for the reference
   - TEXT_FILE: if not empty, output text coverage to this file.
   - STAT_FILE: if not empty, output contig coverage statistics to
this file.
   - HYBRID_FILE: if not empty, output an intermediate file, 
HYBRID_FILE.hybridtemp, for use with HybridSelectionAnalysis.pl.  Outputs 
target window #, # of baits in the target window, and @ reads.  Also, if 
HYBRID_HISTO is true, outputs the file HYBRID_FILE.hybridhistogram that has 
the format:
Window 1
0  1  <# catches base 1>
0  2  <# catches base 2>
...
1  1  <# catches base 1>
1  2  <# catches base 2>
...
The above format is suitable for extraction of specific windows by using
grep.  (This is better than a histogram file for each window as there may be
many of them.)  If one really wants every base printed out, set PRINT_ZEROS
to True.
   - HYBRID_HISTO: see above (default False - this will take up a lot of 
disk space).
   - ZERO_CONTIGS: for coverage stats, print contigs with zero coverage
(default True).
   - HEAD: if not empty, will add num_bases_with_zero_coverage metric 
to the metrics file with prefix HEAD.  SUMMARY must be set to true for
this to apply.  (The value is computed in the summary.)
   - MULTIPLES_RANDOM (default True): if multiple alignment, assign
the read to one position at random . If false, do not lay down the
read.
   - MULTIPLES_ALL (default False): if multiple alignment, assign the
read to all positions meeting MAX_ERRORS.  If false, do not lay down
the read.
   - WIDEN_COVERAGE (default 0): if positive, widen the coverage
     intervals from the alignments by this amount on each side.  If
     negative, narrow them.
   - SUMMARY: print summary coverages to stdout (default True)
   - DELIM: string separating the alignment file names.
   - MAX_ERRORS: do not use alignments with more errors than this. This 
applies even when looking at multiple aligments: we use only the ones that
have <= MAX_ERRORS errors for the random choice.
   - PRINT_ZEROS: print bases with 0 coverage in TEXT_FILE.
   - START_POINT_FILE (default ""): if nonempty, write start point analysis to that file.
   - START_POINT_COUNT_BINS (default "[0,20]"): an IntSet for the bins to use for
   start point pile-height histogram.
   - PROCESS_FW (default True): process alignments that go forward.
   - PROCESS_RC (default True): process alignments that go in reverse.  Both of these options being true means to process all alignments.

If there is more than one alignment file, we expect that each file
corresponds to a separate set of reads and aggregate the coverage from
all files.
*/

#include "MainTools.h"
#include "Basevector.h"
#include "FastaFileset.h"
#include "SeqInterval.h"
#include "CoverageAnalyzer.h"
#include "lookup/LookAlign.h"
#include "lookup/LookAlignFinder.h"
#include "random/Random.h"
#include "TokenizeString.h"
#include "Map.h"
#include "Histogram.h"
#include "solexa/SolexaMetrics.h"

seq_interval Widen(const seq_interval &in, const vec<int> &reflengths, int widen)
{
  seq_interval res = in;
  if (widen!=0) {
    res.SetBegin(min(max(0, res.Begin() - widen), res.End()));
    res.SetEnd(max(min(res.End() + widen, reflengths[in.SeqId()]), res.Begin()));
  }
  return res;
}

typedef StdMap<unsigned int, unsigned int> StartPointCount;
typedef map<unsigned int, unsigned int> PileCount;

// Count the number of occurences of each count in a sequence of
// pair<unsigned int, unsigned int>'s, where the first is the position
// and the second is the number of reads starting there.  That is,
// count the number of times each value of the second occurs, as well
// as the number of 0's between the start positions that occur and
// after the last start position.  Returns the total number of reads,
// which is the sum of all seconds.
template <class In, class Map>
unsigned int countPiles(In curr, In last, Map &counts, unsigned int length)
{
  unsigned int sum=0;
  In next = curr;
  while (curr!=last) {
    // Count the starts pile at *curr
    ++counts[curr->second];
    sum += curr->second;
    ++next;
    // Now count the number of 0's between *curr and *next
    if (next!=last) 
      counts[0] += next->first - curr->first - 1; // still on same target
    else 
      counts[0] += length - curr->first - 1; // changed targets
    curr = next;      // And on to the next
  }
  return sum;
}

int main(int argc, char ** argv) {
  RunTime();
  BeginCommandArguments;
  CommandArgument_StringSet(QLTOUT);
  //CommandArgument_String_OrDefault(REF_LENGTHS, "");
  CommandArgument_String(REF);
  CommandArgument_String_OrDefault_Doc(TEXT_FILE,"","coverage file in one-base-per-line format");
  CommandArgument_String_OrDefault_Doc(STAT_FILE,"","coverage statistics of contigs file in one-contig-per-line format");
  CommandArgument_String_OrDefault_Doc(HYBRID_FILE,"","coverage statistics of windows in hybrid selection");
  CommandArgument_Bool_OrDefault_Doc(HYBRID_HISTO, False, "outputs a histogram of the coverage of the windows");
  CommandArgument_Bool_OrDefault_Doc(ZERO_CONTIGS, True, "print contigs with zero coverage in statistics file");
  CommandArgument_String_OrDefault_Doc(HEAD,"","prefix of the metrics file");
  CommandArgument_String_OrDefault_Doc(COVERED_FILE,"", "coverage summary in HoIntervalWithId format");
  CommandArgument_Bool_OrDefault_Doc(SUMMARY, True, "print summary coverages to stdout");
  CommandArgument_Bool_OrDefault_Doc(MULTIPLES_RANDOM, True,  "if multiple alignment, assign the read to one position at random. If false, do not lay down the read");
  CommandArgument_Bool_OrDefault_Doc(MULTIPLES_ALL, False, "if multiple alignment, assign the read to all positions meeting MAX_ERRORS. If false, do not lay down the read");
  CommandArgument_Int_OrDefault_Doc(MAX_ERRORS, -1, "do not use alignments with more errors than this. This applies even when looking at multiple aligments: we use only the ones that have <= MAX_ERRORS errors for the random choice.");
  CommandArgument_Int_OrDefault_Doc(WIDEN_COVERAGE, 0, "if positive, widen the coverage intervals from the alignments by this amount on each side.  If negative, narrow them.");
  CommandArgument_Bool_OrDefault_Doc(PRINT_ZEROS, False, "print bases with 0 coverage in TEXT_FILE.");
  CommandArgument_String_OrDefault_Doc(START_POINT_FILE,"","write start point analysis to this file");
  CommandArgument_IntSet_OrDefault_Doc(START_POINT_COUNT_BINS, "[0,20]", "an IntSet for the bins to use for start point pile-height histogram");
  CommandArgument_Bool_OrDefault_Doc(PROCESS_FW, True, "Process forward alignments.");
  CommandArgument_Bool_OrDefault_Doc(PROCESS_RC, True, "Process reverse [compliment] alignments.");
  
  EndCommandArguments;

  ForceAssert(!(MULTIPLES_RANDOM && MULTIPLES_ALL)); // These can both be false, but both true is meaningless

  ForceAssert(PROCESS_FW || PROCESS_RC);  // One of them has to be true, if not both.

  vec<int> * refLengths = new vec<int>;//CoverageAnalyzer takes ownership
  if (REF.Contains(".fastb")) {
    vecbasevector ref;
    ref.ReadAll(REF);
    transform(ref.begin(), ref.end(), back_inserter(*refLengths),
	      mem_fun_ref(&basevector::size));
  } 
  else { //do not read in the whole file to avoid memory overload.
    temp_file fname("reflengths_XXXXXXX");
    SystemSucceed("FastaSizes < " + REF + " > " + fname);
    {
      Ifstream(is, fname);
      int size;
      while (true) { is >> size; if (!is) break; refLengths->push_back(size); }
    }
  }
  const int nRefs = refLengths->size();
  longlong refSize = 0;
  for (int i=0; i<nRefs; ++i)
    refSize += (*refLengths)[i];

  vec<String> alignfiles(QLTOUT);

  vec<seq_interval> intervals;
  order_lookalign_QueryErrors sorter;
  vec<StartPointCount> startCounts_fw(nRefs), startCounts_rc(nRefs);
  // Below needed for keeping track of the number of catches 
  // per window (reference) used in hybrid selection
  vector<unsigned int> countCatches(nRefs, 0);
  const int dotter = 10000;
  int firstchoice, lastchoice;
  int nReads = 0;
  for (int i=0; i != alignfiles.isize(); ++i) {
    cout << "Loading data from file " << alignfiles[i] << endl;
    cout << "Processing aligned reads in groups of  " << dotter << endl;
    int dotcount=0;
    for (LookAlignFinder finder(alignfiles[i]); !finder.empty(); ++finder) {
      DotMod(cout, ++dotcount, dotter);
      int r = finder.QueryId();
      vec<look_align_plus> & aligns = finder.Aligns();
      sort(aligns.begin(), aligns.end(), sorter);
      firstchoice = 0;
      lastchoice = 0;
      switch (aligns.size()) {
	// Determine our choice(s) of alignment
      case 0:
	break;
      case 1:
	lastchoice = 1;
	break;
      default: // More than one plausible alignment.
	if (MULTIPLES_RANDOM) { // pick a random alignment of the read...
	  int maxchoice=-1;
	  if (MAX_ERRORS < 0) maxchoice = aligns.size();
	  else {
	    for (int m = aligns.size()-1; m >= 0; --m) {
	      if  (aligns[m].Errors() <= MAX_ERRORS) {
		maxchoice = m+1;
		break;
	      }
	    }
	  }
	  if (maxchoice >= 0) {
	    firstchoice = randomx() % maxchoice;
	    lastchoice = 1+firstchoice;
	  }
	} else if (MULTIPLES_ALL) {
	  lastchoice = aligns.size();
	}
	break;
      }
      for (int choice=firstchoice; choice<lastchoice; ++choice) {
	// Process the choice(s), uniformly
	const look_align_plus & la = aligns[choice];
	if (MAX_ERRORS<0 || la.Errors() <= MAX_ERRORS) {
         if ((!la.rc1 && PROCESS_FW) || (la.rc1 && PROCESS_RC)) {
	  // Here's what we do for each chosen accepted alignment:
	  ++nReads;
	  ++countCatches[la.target_id];
	  intervals.push_back(Widen(la, *refLengths, WIDEN_COVERAGE));
	  if (!START_POINT_FILE.empty()) {
	    if (la.rc1) {
	      ++startCounts_rc[la.target_id][la.Pos2()-1];
	    } else {
	      ++startCounts_fw[la.target_id][la.pos2()];
	    }
	  }
         }
	}
      }
    }
  }
  cout << "\n"; //clean up after dots

  PRINT(nReads);
  
  CoverageAnalyzer cov(intervals, refLengths);

  if (SUMMARY) {
    vec<longlong> counts;
    cov.CountAllCoverages(counts);
    cout << "Bases without coverage: " << counts[0] << endl;
    cout << "\nCoverage\t#bases\n";
    for (int i=0; i != counts.isize(); ++i) {
      if (counts[i]) {
	cout << i << "\t" << counts[i] << "\n";
      }
    }
    if (!HEAD.empty()) {
      solexa_metric_db db(HEAD + ".metrics");
      db.SetValue("num_bases_with_zero_coverage", counts[0]);
      db.WriteMetrics(HEAD + ".metrics");
    }
  }

  if (!TEXT_FILE.empty()) {
    Ofstream(text, TEXT_FILE);
    cov.PrintTextCoverage(text, true, PRINT_ZEROS);
  }

  if (!STAT_FILE.empty()) {
    Ofstream(statfile, STAT_FILE);
    cov.PrintCoverageStats(statfile, ZERO_CONTIGS);
  }

  if (!HYBRID_FILE.empty()) {
    Ofstream(hybridfile, HYBRID_FILE + ".hybridtemp");
    int catchSize = countCatches.size();
    for (int i = 0; i < catchSize; ++i)
       hybridfile << i << "\t" << countCatches[i] << endl; 
    if (HYBRID_HISTO) {
       Ofstream(hybridhisto, HYBRID_FILE + ".hybridhistogram");
       cov.PrintTextCoverage(hybridhisto, true, PRINT_ZEROS);
    }
  }

  if (!COVERED_FILE.empty()) {
    vec<seq_interval> covered;
    cov.GetCoveragesAtLeast(1, covered);
    Ofstream(file, COVERED_FILE);
    HoIntervalWithId ho;
    for (unsigned int i=0; i<covered.size(); ++i){
      ho.id = covered[i].SeqId();
      ho.Set(covered[i].Begin(), covered[i].End());
      file << ho << "\n";
    }
  }

  if (!START_POINT_FILE.empty()) {
    const double readRate = double(nReads)/double(refSize);
    PRINT(readRate);
    Ofstream(out, START_POINT_FILE);
    histogram<int> hist;
    hist.AddBins(START_POINT_COUNT_BINS);
    out << "Contig\t#bases\t#reads\treadspb\trdcov\t#fw reads\t#rc reads\tz_fwrc\tdup\txdup"
	<< "\tSPHist";
    for (int i=0; i<START_POINT_COUNT_BINS.isize(); ++i) {
      out << "\t" << i;
    }
    out << "\n";
    for (int c=0; c<refLengths->isize(); ++c) {
      PileCount piles;
      hist.ResetCounts();
      const int R = (*refLengths)[c];
      int fw = countPiles(startCounts_fw[c].begin(), startCounts_fw[c].end(), piles, R);
      int rc = countPiles(startCounts_rc[c].begin(), startCounts_rc[c].end(), piles, R);

      // Record in histogram the number of piles of each size
      int dups = 0; 
      for (PileCount::iterator it=piles.begin(); it!= piles.end(); ++it) {
	hist.AddIdenticalData(it->first, it->second);
	if (it->first > 1)
	  dups += (it->first - 1) * (it->second);
      }

      // Compute some basic statistics
      const double N=fw+rc, M = max(fw,rc), arg = (2*M-N)/sqrt(N);
      const double rate = N/R;

      // Compute expectation of #duplicates.  The effective reference
      // size is 2R because fw and rc reads are treated separately.
      int i=0;
      double p=Poisson(rate/2.0, i), expdups = 0;
      while (p>0) {
	if (i>1) expdups += 2*R*p*(i-1);
	++i;
	p=Poisson(rate/2.0, i);
      }
      
      const char goodRate = (rate>5*readRate ? 'H' : (rate < readRate/5 ? 'L' : '-'));
      out << c << "\t" << R << "\t" << N << "\t" << rate << "\t" << goodRate << "\t"
	  << fw << "\t" << rc << "\t" << arg << "\t" << dups << "\t" << (dups-expdups) << "\t\t";
      hist.PrintDataAsLine(out);
    }
  }

  return 0;
}



  
