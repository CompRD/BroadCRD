// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef PREDIGEST_ALIGNMENTS_UTIL_H
#define PREDIGEST_ALIGNMENTS_UTIL_H

#include "String.h"
#include "lookup/LookAlign.h"



/*
 * TargetId
 *
 * Find the target_id of a look_align (which is passed as a String).
 */
int TargetId( const String &lookalign )
{
  return( lookalign.After( "\t" ).After( "\t" ).After( "\t" ).After( "\t" )
	  .After( "\t" ).After( "\t" ).Before( "\t" ).Int( ) );
}



/*
 * SameQueryId
 * ordering functor
 *
 * Check if two look_align objects have the same query_id. Notice that
 * they are passed as Stringes.
 */
struct SameQueryId
  : public binary_function<const String&, const String&, bool>
{
public:
  bool operator() ( const String &left, const String &right ) {
    int left_id = left.After( "\t" ).Before( "\t" ).Int( );
    int right_id = right.After( "\t" ).Before( "\t" ).Int( );
    
    return ( left_id == right_id );
  }
};



/*
 * SortByQueryId
 * ordering functor
 *
 * Sort two look_align objects (as Stringes) by their query id. 
 */
struct SortByQueryId
  : public binary_function<const String&, const String&, bool>
{
public:
  bool operator() ( const String &left, const String &right ) {
    int left_id = left.After( "\t" ).Before( "\t" ).Int( );
    int right_id = right.After( "\t" ).Before( "\t" ).Int( );
    
    return ( left_id < right_id );
  }
};



/*
 * SortByStartOnTarget
 * ordering functor
 *
 * Sort two look_align objects by their starting point on the target_id.
 * Look_aligns are passed as Stringes.
 */
struct SortByStartOnTarget
  : public binary_function<const String&, const String&, bool>
{
public:
  bool operator() ( const String &left, const String &right ) {
    int left_id = left.After( "\t" ).After( "\t" ).After( "\t" )
      .After( "\t" ).After( "\t" ).After( "\t" ).After( "\t" )
      .Before( "\t" ).Int( );
    
    int right_id = right.After( "\t" ).After( "\t" ).After( "\t" )
      .After( "\t" ).After( "\t" ).After( "\t" ).After( "\t" )
      .Before( "\t" ).Int( );
    
    return ( left_id < right_id );
  }
};



#endif
