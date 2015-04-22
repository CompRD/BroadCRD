// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
//

// Consensus-vs-Known Differences basic statistics.

#ifndef CK_DIFF_H
#define CK_DIFF_H

#include "system/System.h"



/*
 * class ck_diff
 *
 * This data-stracuture manages basic statistics for differences between
 * a contig (consensus) and a given known sequence.
 */
class ck_diff {

public:
  
  ck_diff ( );
  
  int TotalDiff( ) const;

  void PrintCompact( ostream &out, bool oneline = true ) const;

  void Print( ostream &out, bool oneline = false, bool newline = true ) const;
  
  ck_diff &operator+= ( const ck_diff &addendum );
  
  
public:

  int n_total;         // total number of bases analyzed
  int n_pads_k;        // number of indels as pads_on_known
  int n_pads_c;        // number of indels as pads_on_consensus
  int n_confirm_c;     // evidence to support consensus only (count)
  int n_confirm_k;     // evidence to support known only
  int n_confirm_both;  // evidence to support both known and consensus
  int n_confirm_none;  // neither known nor consensus are supported

};



#endif
