/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Alignment.h"
#include "Basevector.h"

// This class generates all the perfect alignments between sequences
// in some vecbasevector.

class PerfectAlignerImp;

class PerfectAligner {
 public:
  enum Behavior {
    findProperOnly,
    findSemiproper,
    findImproper
  };

  PerfectAligner( int kmerSize, 
                  Behavior behavior,
                  ostream* pLog = 0 );
  ~PerfectAligner();

  void SetKmerSize( const int kmerSize );
  void SetBehavior( const Behavior behavior );
  void SetMaxKmerFreq( const int maxKmerFreq );

  // Calculate all the perfect alignments conforming to the current Behavior.
  // If partition is less than zero, compare all sequences to all sequences.  
  // Otherwise, compare all the sequences with ids less than partition to
  // all the sequences with ids greater than or equal to partition.
  void Align( const vecbasevector& sequences, 
              vec<alignment_plus>& perfectAligns,
              const int partition = -1,
	      const bool append = false );

 private:
  int m_kmerSize;
  int m_maxKmerFreq;
  Behavior m_behavior;
  ostream* m_pLog;

  PerfectAlignerImp* m_pImp;
};

