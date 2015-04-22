// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef BUILD34_H
#define BUILD34_H

#include "system/Assert.h"
#include "String.h"



/*
 * ToStandard
 *
 * Converts build34 entry positions to standard chromosome names. Assert
 * if pos is not a valid entry.
 */
String ToStandard( int pos )
{
  if ( pos == 0 )
    return "DR51";
  else if ( 0 < pos && pos < 23 )
    return ToString( pos );
  else if ( 23 == pos)
    return "X";
  else if ( 24 == pos )
    return "Y";
  else if ( 25 == pos  )
    return ( "_unplaced" );
  
  ForceAssert( 1 == 0 );
  return ( "" );
}



/*
 * ToBuild34
 *
 * Converts standard chromosome names to build34 entry positions. Assert if
 * chr_name is not in the valid form.
 */
int ToBuild34( const String &chr_name )
{
  if ( chr_name.IsInt( ) ) {
    int pos = chr_name.Int( );
    
    if ( 0 < pos && pos < 23 )
      return pos;
    else
      ForceAssert( 1 == 0 );
  }
  else if ( chr_name == "X" ) {
    return 23;
  }
  else if ( chr_name == "Y" ) {
    return 24;
  }
  else if ( chr_name.Contains( "_unplaced", 0 ) ) {
    return 25;
  }
  else if ( chr_name.Contains( "DR51", 0 ) ) {
    return 0;
  }    
  
  ForceAssert( 1 == 0 );
  return -1;
}



#endif
