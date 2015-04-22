// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef PROTO_CONTIG_ID_H
#define PROTO_CONTIG_ID_H

#include <map>

#include "Vec.h"



/*
 * class proto_contig_id
 *
 * A container for a contig identifier. To be used with the OnMaster code
 * (chimp on humana, etc.)
 */
class proto_contig_id {
  
public:

  proto_contig_id( );
  
  proto_contig_id( int cg_id, int chr_id, int pos );
  
  int ContigId( ) const;
  
  int ChrId( ) const;
  
  int ChrPos( ) const;
  
  friend bool operator==( const proto_contig_id &l, const proto_contig_id &r );

  friend istream &operator>> ( istream &in, proto_contig_id &id );

  friend ostream &operator<< ( ostream &out, const proto_contig_id &outx );
  
  
private:

  int contig_id_;  // contig id
  int chr_id_;     // chromosome id
  int chr_pos_;    // position in chromosome

};



/*
 * order_protocontigid_Contig
 * ordering functor
 */
struct order_protocontigid_Contig :
  public binary_function<const proto_contig_id&, const proto_contig_id&, bool>
{
public:
  bool operator( ) ( const proto_contig_id &l, const proto_contig_id &r ) {
    return ( l.ContigId( ) < r.ContigId( ) );
  }
};



/*
 * order_protocontigid_ChrPos
 * ordering functor
 */
struct order_protocontigid_ChrPos :
  public binary_function<const proto_contig_id&, const proto_contig_id&, bool>
{
public:
  bool operator( ) ( const proto_contig_id &l, const proto_contig_id &r ) {
    if ( l.ChrId( ) == r.ChrId( ) )
      return ( l.ChrPos( ) < r.ChrPos( ) );
    return ( l.ChrId( ) < r.ChrId( ) );
  }
};



#endif
