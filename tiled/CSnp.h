/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef C_SNP_H
#define C_SNP_H

#include "String.h"
#include "tiled/CharBaseUtil.h"

/**
 * class CSnp
 *
 * Container for bases/scores for a SNP event. Types are defined in
 * tiled/CharBaseUtil (pads are allowed bases), and quality scores are
 * capped.
 */
class CSnp {
  
public:

  CSnp( );
  
  CSnp( char b1, int q1, char b2 = empty_base, int q2 = empty_qual );
  
  void Set( char b1, int q1, char b2 = empty_base, int q2 = empty_qual );
  
  // Unpack (snp_base,snp_qual) onto (b1_,q1_) and (b2_,q2_).
  //  void Unpack( char snp_base, int snp_qual );
  
  // Pack (b1_,q1_) and (b2_,q2_) into (snp_base,snp_qual).
  void Pack( char &snp_base, int &snp_qual, const int *cap = 0 ) const;

  // Empty SNP (uninitialized).
  bool IsEmptySnp( ) const;

  // Find (b1_,q1_). If second=true, find (b2_,q2_).
  void Base( char &base, int &qual, bool second = false ) const;
  
  // Print (b1_,q1_). If second=true print (b2_,q2_).
  String PrintBase( bool second = false ) const;
  
  // Print full snp.
  String PrintSnp( ) const;

  // Two CSnps are equal iff the bases match (quals may differ).
  friend bool operator== ( const CSnp &left, const CSnp &right );

  // Negation of the above (for convenience).
  friend bool operator!= ( const CSnp &left, const CSnp &right );

  // This will assert if left and right are not equal.
  friend CSnp operator+ ( const CSnp &left, const CSnp &right );
  
  
private:
  
  // Cap quals, upper case all bases. Store in (b1,q1) the highest qual.
  void Standardize( );
  
  
private:
  
  char b1_;  // first base
  int q1_;   // first qual
  char b2_;  // second base (may be empty_base)
  int q2_;   // second qual (may be empty qual)
  
};

#endif
