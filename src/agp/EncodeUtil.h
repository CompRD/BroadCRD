// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
//

// encode_hit: data structures encapsulating agp-related info.

#ifndef ENCODE_UTIL_H
#define ENCODE_UTIL_H

#include <map>

#include "CoreTools.h"
#include "FastIfstream.h"
#include "math/Functions.h"
#include "Qualvector.h"
#include "Superb.h"
#include "TokenizeString.h"
#include "agp/AgpFile.h"



/*
 * class encode_hit
 *
 * An encode region is identified by chromosome name, and by a region given
 * as a window on the chromosome.
 */
class encode_hit {
  
public:

  encode_hit( ) { }
  
  // Format for desname: name:start-stop (e.g. "chrUn:5000-5999").
  encode_hit( const String &desname ) {
    vec<char> seps;
    seps.push_back( '-' );
    seps.push_back( ':' );
    
    vec<String> tokens;
    Tokenize( desname, seps, tokens );
    ForceAssert( tokens.size( ) == 3 );

    chr_name_ = tokens[0];
    begin_ = tokens[1].Int( );
    end_ = 1 + tokens[2].Int( );
  }

  String ChrName( ) const { return chr_name_; }

  int Begin( ) const { return begin_; }
  
  int End( ) const { return end_; }
  
  String PrettyName( ) const {
    String str_start = ToString( begin_ );
    String str_stop = ToString( end_ - 1 );
    return chr_name_ + ":" + str_start + "-" + str_stop;
  }

  friend bool operator< ( const encode_hit &left, const encode_hit &right ) {
    if ( left.chr_name_ == right.chr_name_ )
      return ( left.begin_ < right.begin_ );
    return ( left.chr_name_ < right.chr_name_ );
  }


private:

  String chr_name_;
  int begin_;
  int end_;
  
};



/*
 * Interesects
 *
 * Return the id of the encode_hit intersecting the region determined by
 * chr_name, begin, and end. Return -1 if no such region is found.
 */
int Intersects( const vec<encode_hit> &hits,
		const String &chr_name,
		int begin,
		int end )
{
  for (int hit_id=0; hit_id<(int)hits.size( ); hit_id++) {
    if ( hits[hit_id].ChrName( ) != chr_name )
      continue;
    int left_int = Max( begin, hits[hit_id].Begin( ) );
    int right_int = Min( end, hits[hit_id].End( ) );
    if ( right_int - left_int > 0 )
      return hit_id;
  }
  return -1;
}



#endif
