// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#include "String.h"
#include "TokenizeString.h"
#include "tiled/PrimateContigLogParser.h"



/*
 * pclog_parser
 * Constructor
 */
pclog_parser::pclog_parser ( ) :
  chr_id_ ( -1 ),
  contig_id_ ( -1 ),
  n_reads_ ( -1 ),
  begin_ ( -1 ),
  len_ ( -1 )
{ }



/*
 * pclog_parser
 * Constructor
 */
pclog_parser::pclog_parser( int chr_id,
			    int contig_id,
			    int n_reads,
			    int begin,
			    int len ) :
  chr_id_ ( chr_id ),
  contig_id_ ( contig_id ),
  n_reads_ ( n_reads ),
  begin_ ( begin ),
  len_ ( len )
{ }



/*
 * pclog_parser
 * Set
 */
void pclog_parser::Set( int chr_id,
			int contig_id,
			int n_reads,
			int begin,
			int len )
{
  chr_id_ = chr_id;
  contig_id_ = contig_id;
  n_reads_ = n_reads;
  begin_ = begin;
  len_ = len;
}



/*
 * pclog_parser
 * operator<<
 */
ostream &operator<< ( ostream &out, const pclog_parser &pclog )
{
  out << pclog.chr_id_ << "\t"
      << pclog.contig_id_ << "\t"
      << pclog.n_reads_ << "\t"
      << pclog.begin_ << "\t"
      << pclog.len_;

  return out;
}



/*
 * pclog_parser
 * operator>>
 */
istream &operator>> ( istream &in, pclog_parser &pclog )
{
  in >> pclog.chr_id_
     >> pclog.contig_id_
     >> pclog.n_reads_
     >> pclog.begin_
     >> pclog.len_;

  return in;
}



