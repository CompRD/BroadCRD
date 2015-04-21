// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef TILES_CONTIG_PARSER_H
#define TILES_CONTIG_PARSER_H

#include "Basevector.h"
#include "Qualvector.h"
#include "String.h"
#include "Vec.h"
#include "tiled/PaddedSeq.h"



/*
 * class tc_parser
 *
 * Parse ReadsTiling type consensus bases and quality scores. It reads a
 * pair of files at a time (bases and quality scores), appending bases
 * and quals to member vecbasevector and vecqualvector.
 */
class tc_parser {
  
public:
  
  tc_parser( ) { }
  
  void Reserve( int n_reads, longlong n_bases );
  
  void Append( const String &in_fasta, const String &in_qual );
  
  void Write( const String &out_fastb,
	      const String &out_qualb,
	      const String &out_pads ) const;
	      
	      
  
  
private:
  
  vecbasevector bases_;
  vecqualvector quals_;
  vec<padded_seq> pads_;
  
};



#endif
