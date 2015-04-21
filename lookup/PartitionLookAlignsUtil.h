// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef PARTITION_LOOK_ALIGNS_UTIL_H
#define PARTITION_LOOK_ALIGNS_UTIL_H

#include "String.h"



/*
 * order_Filename
 *
 * This is a very local ordering functor, which works only on file names
 * of the type "xxx.yyy.aaa.bbb", where aaa < bbb are two positive
 * integers. It sorts Stringes by the aaa's.
 */
struct order_Filename
  : public binary_function<const String&, const String&, bool>
{
  bool operator() ( const String &left, const String &right ) {
    int n_left = left.After( "." ).After( "." ).Before( "." ).Int( );
    int n_right = right.After( "." ).After( "." ).Before( "." ).Int( );
    return ( n_left < n_right );
  }
};



/*
 * order_QueryId
 * ordering functor
 *
 * Sort two look_align objects (as Stringes) by their query id. 
 */
struct order_QueryId
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
 * GetRange
 *
 * It returns the range [nnn, mmm) in a file of the type xxx.yyy.nnn.mmm
 * as explained above.
 */
pair<int, int> GetRange( const String& file_name )
{
  int left = file_name.After( "." ).After( "." ).Before( "." ).Int( );
  int right = file_name.After( "." ).After( "." ).After( "." ).Int( );
  
  return ( make_pair( left, right ) );
}



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



#endif
