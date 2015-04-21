/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef CORRECTOR_H
#define CORRECTOR_H

#include "String.h"
#include "tiled/CharBaseUtil.h"
#include <istream>
#include <ostream>

/**
 * class corrector
 *
 * A simple structure encoding what and how to correct in a given
 * genomic sequence. See tiled/CharBaseUtil.h for legal values for
 * bases and quals (represented as char's and int's).
 */
class corrector {
  
public:
  
  corrector( );
  
  corrector( int id, int pos,
	     char oldB, char newB,
	     int oldQ, int newQ );

  corrector( int id, int pos,
	     char oldB, char newB_0, char newB_1,
	     int oldQ, int newQ_0, int newQ_1 );

  String OldBase( ) const;
  String NewBase( ) const;
  char CharOldBase( ) const { return old_base_; }
  char CharNewBase( ) const { return new_base0_; }
  char CharNewBase1( ) const { return new_base1_; }
  
  void BinaryOut( std::ostream &out ) const;
  
  void BinaryIn( std::istream &in );
  
  friend std::ostream &operator<< ( std::ostream &out, const corrector &crc );

  friend std::istream &operator>> ( std::istream &in, corrector &crc );
  
  
public:
  
  int id_;            // sequence identifier
  int pos_;           // position in sequence
  char old_base_;     // old base
  char new_base0_;    // new base
  char new_base1_;    // optional second haplotype
  int old_qual_;      // old qual as an int
  int new_qual0_;     // new qual as an int
  int new_qual1_;     // optional second haplotype qual as an int
  
};

#endif
