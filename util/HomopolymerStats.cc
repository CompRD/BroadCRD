// Copyright (c) 2004 Broad Institute of MIT and Harvard
//
// HomopolymerStats.cc
//

// This program takes as input a fastb file of bases,
// and computes the percentage of bases which fall in
// homopolymeric stretches of lengths 1,2,3,4,5,6 and >=7.


// Command line options:
//   FASTB_FILE
//   INTERVALS (optional), HoIntervalWithID format.

#include "MainTools.h"
#include "Vec.h"
#include "Basevector.h"
#include "math/HoInterval.h"
#include "FastaFileset.h"

//----------------------------------------------------------------------

int main( int argc, char *argv[] )
{
  //------------------------------
  // Initialzation
  //

  RunTime();

  BeginCommandArguments;

  // input file of basevectors whose homopolymer statistics are to be computed
  CommandArgument_String_Abbr(FILE, F);
  // Optional HoIntervalWithId-format vector of intervals to restrict processing to.
  CommandArgument_String_OrDefault(INTERVALS, "");

  EndCommandArguments;


  //------------------------------
  // Input data
  //

  cout << Date() << ": Reading " << FILE << "..." << endl;
  vecbasevector  contigs;
  LoadReads(contigs, FILE);
  cout << Date() << ": Reading " << FILE << "...done." << endl;

  
  vec<HoIntervalWithId> intervals;
  if (!INTERVALS.empty()) {
    Ifstream(in, INTERVALS);
    while (in) {
      HoIntervalWithId tmp;
      in >> tmp;
      if (in)
	intervals.push_back(tmp);
    }
  } else {
    // Create intervals matching contigs so following code doesn't
    // have to be duplicated
    for (size_t i=0; i<contigs.size(); ++i)
      intervals.push_back(HoIntervalWithId(0, contigs[i].isize(), i));
  }

  //------------------------------
  // Compute homopolymer statistics
  //

  const int N = intervals.size();
  cout << Date() << ": Computing homopolymer statistics for "
       << N << " intervals..." << endl;

  const int maxMult = 8;
  vec<longlong> freq_of_mult(maxMult,0);
  longlong num_bases = 0;

  for (int i = 0; i < N; i++) {
    const basevector & contig = contigs[intervals[i].id];
    const int end = intervals[i].Stop();
    ForceAssertLe(end, contig.isize());
    num_bases += intervals[i].Length();

    int pos=intervals[i].Start();
    while (pos < end) {
      const unsigned char curr_base = contig[pos];
      int mult = 0;

      do {
	++pos;
	++mult;
      } while (pos < end && contig[pos] == curr_base);

      freq_of_mult[(mult<maxMult ? mult : maxMult-1)] += mult;
    }
  }

  cout << Date() << ": done." << endl;

  //------------------------------
  // Output results
  //

  cout << endl;

  cout.setf(ios::fixed);
  cout.setf(ios::right);
  cout.setf(ios::showpoint);


  cout << "mult  frequency" << endl;
  for (int i = 1; i < maxMult-1; i++)
    {
      cout << "  " << setw(1) << i
	   << "  " << setw(10) << freq_of_mult[i]
	   << " (" << setw(4) << setprecision(1)
	   << static_cast<double>(freq_of_mult[i] * 100) / static_cast<double>(num_bases)
	   << "%)" << endl;
    }
  const int i = maxMult-1;
  cout << ">=" << setw(1) << i
       << "  " << setw(10) << freq_of_mult[i]
       << " (" << setw(4) << setprecision(1)
       << static_cast<double>(freq_of_mult[i] * 100) / static_cast<double>(num_bases)
       << "%)" << endl;

} // main()



