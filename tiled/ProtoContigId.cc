// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#include "Vec.h"
#include "tiled/ProtoContigId.h"



/*
 * proto_contig_id
 * Constructor
 */
proto_contig_id::proto_contig_id( )
{ }



/*
 * proto_contig_id
 * Constructor
 */
proto_contig_id::proto_contig_id( int cg_id, int chr_id, int pos ) :
  contig_id_ ( cg_id ),
  chr_id_ ( chr_id ),
  chr_pos_ ( pos )
{ }



/*
 * proto_contig_id
 * ContigId
 */
int proto_contig_id::ContigId( ) const
{
  return contig_id_;
}



/*
 * proto_contig_id
 * ChrId
 */
int proto_contig_id::ChrId( ) const
{
  return chr_id_; 
}



/*
 * proto_contig_id
 * ChrPos
 */
int proto_contig_id::ChrPos( ) const
{
  return chr_pos_;
}



/*
 * proto_contig_id
 * operator==
 */
bool operator== ( const proto_contig_id &left, const proto_contig_id &right )
{
  return ( left.contig_id_ == right.contig_id_ &&
	   left.chr_id_ == right.chr_id_ &&
	   left.chr_pos_ == right.chr_pos_ );
}



/*
 * proto_contig_id
 * operator>>
 */
istream &operator>> ( istream &in, proto_contig_id &id )
{
  in >> id.contig_id_
     >> id.chr_id_
     >> id.chr_pos_;
  
  return in;
}



/*
 * proto_contig_id
 * operator<<
 */
ostream &operator<< ( ostream &out, const proto_contig_id &id )
{
  out << id.contig_id_ << "\t"
      << id.chr_id_ << "\t"
      << id.chr_pos_ << "\n";
  
  return out;
}



