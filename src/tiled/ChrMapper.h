// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef CHR_MAPPER_H
#define CHR_MAPPER_H

#include <iostream>

using namespace std;



/*
 * class chr_mapper
 *
 * Data class for the position of a contig on a chromosome.
 */
class chr_mapper {
  
public:
  
  chr_mapper( );

  chr_mapper( int chr, int id, int pos, int begin );
  
  void Set( int chr, int id, int pos, int begin );
  
  void SetChr( int chr ) { chr_ = chr; }

  void SetId( int id ) { id_ = id; }

  void SetPos( int pos ) { pos_ = pos; }

  void SetBegin( int begin ) { begin_ = begin; }

  int Chr( ) const { return chr_; }
  
  int Id( ) const { return id_; }
  
  int Pos( ) const { return pos_; }
  
  int Begin( ) const { return begin_; }

  friend bool operator< ( const chr_mapper &left, const chr_mapper &right );

  friend istream &operator>> ( istream &in, chr_mapper &cmap );
  
  friend ostream &operator<< ( ostream &out, const chr_mapper &cmap );
  
  
private:
  
  int chr_;     // chromosome id
  int id_;      // contig id
  int pos_;     // position on chromosome
  int begin_;   // begin on chromosome
  
};



#endif
