/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__ASSISTED__C_LIB_LINKS__H
#define PATHS__ASSISTED__C_LIB_LINKS__H

#include "Vec.h"

/**
 * class CRawLink
 *
 * One separation between two oriented contigs. It uses even ids for
 * fw contigs, odd ids for rc contigs.
 */
class CRawLink {
  
public:

  CRawLink( ) :
    id1_ ( -1 ), id2_ ( -1 ), sep_ ( 0 ) { }

  CRawLink( int id1, int id2, int sep ) :
    id1_ ( id1 ), id2_ ( id2 ), sep_ ( sep ) { }
  
  int Id1( ) const { return id1_; }

  int Id2( ) const { return id2_; }

  int Sep( ) const { return sep_; }
  
  pair<int,int> Ids( ) const { return make_pair( id1_, id2_ ); }

  friend ostream &operator<< ( ostream &out, const CRawLink &link ) {
    out << link.id1_ / 2 << ( link.id1_ % 2 == 0 ? "+" : "-" ) << "\t"
	<< link.id2_ / 2 << ( link.id2_ % 2 == 0 ? "+" : "-" ) << "\t"
	<< link.sep_;
    return out;
  }

  friend bool operator< ( const CRawLink &left, const CRawLink &right ) {
    if ( left.id1_ < right.id1_ ) return true;
    if ( left.id1_ > right.id1_ ) return false;
    if ( left.id2_ < right.id2_ ) return true;
    if ( left.id2_ > right.id2_ ) return false;
    return ( left.sep_ < right.sep_ );
  }
  
  
private:

  int id1_;
  int id2_;
  int sep_;

};

/**
 * class CLibLinks
 *
 * All links between two oriented contigs from one single library. Dev
 * is set to 0 if the gaps are computed from the alignments of contigs
 * onto a reference (this is treated as a special library).
 */
class CLibLinks {

public:

  CLibLinks( ) : sep_ ( 0 ), dev_ ( 0 ), name_ ( "" ) { }

  void SetLib( int sep, int dev, const String &name ) {
    sep_ = sep;
    dev_ = dev;
    name_ = name;
  }
  
  void Clear( ) { links_.clear( ); }

  void Reserve( int n_links ) { links_.reserve( n_links ); }

  void AddLink( const CRawLink &link ) { links_.push_back( link ); }

  void SortLinks( ) { sort( links_.begin( ), links_.end( ) ); }

  void BuildFirst( int ntigs ) {
    ForceAssert( this->IsSorted( ) );
    first_.resize( 2 * ntigs );
    for (int ii=links_.isize( )-1; ii>=0; ii--)
      first_[links_[ii].Id1( )] = ii;
  }

  const CRawLink &operator[] ( size_t ii ) const { return links_[ii]; }

  int LibSep( ) const { return sep_; }

  int LibDev( ) const { return dev_; }

  bool IsSorted( ) const { return is_sorted( links_.begin( ), links_.end( ) ); }

  int First( size_t ii ) const {
    ForceAssert( first_.size( ) > 0 );
    return first_[ii];
  }

  vec<int> AllSeps( const int v, const int w ) const {
    vec<int> seps;

    pair<int,int> range = this->RangeOf( v, w );
    int n_links = range.second - range.first;
    if ( range.first == range.second ) return seps;
    
    seps.reserve( n_links );
    for (int ii=range.first; ii<range.second; ii++)
      seps.push_back( links_[ii].Sep( ) );

    return seps;
  }
  
  void AllMates( vec< pair<int,int> > &mates ) const {
    mates.clear( );
    if ( links_.size( ) < 1 ) return;
    int nmates = 1;
    for (size_t ii=1; ii<links_.size( ); ii++)
      if ( links_[ii].Ids( ) != links_[ii-1].Ids( ) )
	nmates++;
    mates.reserve( nmates );
    mates.push_back( links_[0].Ids( ) );
    for (size_t ii=1; ii<links_.size( ); ii++)
      if ( links_[ii].Ids( ) != links_[ii-1].Ids( ) )
	mates.push_back( links_[ii].Ids( ) );
  }
  
  pair<int,int> RangeOf( int v, int w ) const {
    int beg = first_[v];
    if ( beg < 0 )
      return make_pair( 0, 0 );
    
    while ( links_[beg].Id1( ) == v && links_[beg].Id2( ) != w ) beg++;
    if ( links_[beg].Id1( ) != v || links_[beg].Id2( ) != w )
      return make_pair( 0, 0 );
    
    int end = beg + 1;
    while( links_[end].Id1( ) == v && links_[end].Id2( ) == w ) end++;
    if ( links_[end-1].Id1( ) != v || links_[end-1].Id2( ) != w )
      return make_pair( 0, 0 );

    return make_pair( beg, end );
  }

  void PrintLinks( ostream &out ) const {
    ForceAssert( is_sorted( links_.begin( ), links_.end( ) ) );
    out << "IMPLIED GAPS FOR LIB " << name_ << " (dev: " << dev_ << ")\n\n";
    for (size_t ii=0; ii<links_.size( ); ii++)
      out << links_[ii] << "\n";
    out << endl;
  }
  
  
private:

  int sep_;              // library given mean sep
  int dev_;              // library given mean dev
  String name_;          // name of library
  vec<CRawLink> links_;  // all separations
  vec<int> first_;       // first_[n]: first link with id1 = n (or -1)

};

#endif
