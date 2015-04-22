// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef BUILD33_H
#define BUILD33_H

#include "system/Assert.h"
#include "String.h"



/*
 * ToStandard
 *
 * Converts build33 entry positions to standard chromosome names. Assert
 * if pos is not a valid entry.
 */
String ToStandard( int pos )
{
  if ( pos == 0 )
    return "1";
  else if ( 0 < pos && pos < 11 )
    return ToString( pos + 9 );
  else if ( pos == 11 )
    return "2";
  else if ( 11 < pos && pos < 15 )
    return ToString( pos + 8 );
  else if ( 14 < pos && pos < 22 )
    return ToString( pos - 12 );
  else if ( 22 == pos )
    return "X";
  else if ( 23 == pos )
    return "Y";
  else if ( pos < 162 )
    return ( "_unplaced" + ToString( pos - 24 ) );

  ForceAssert( 1 == 0 );
  return ( "" );
}



/*
 * ToBuild33
 *
 * Converts standard chromosome names to build33 entry positions. Assert if
 * chr_name is not in the valid form.
 */
int ToBuild33( const String &chr_name )
{
  if ( chr_name.IsInt( ) ) {
    int pos = chr_name.Int( );
    
    if ( pos == 1 )
      return 0;
    if ( 9 < pos && pos < 20 )
      return pos - 9;
    else if ( pos == 2 )
      return 11;
    else if ( 19 < pos && pos < 23 )
      return pos - 8;
    else if ( 2 < pos && pos < 10 )
      return pos + 12;
    else
      ForceAssert( 1 == 0 );
  }
  else if ( chr_name == "X" ) {
    return 22;
  }
  else if ( chr_name == "Y" ) {
    return 23;
  }
  else if ( chr_name.Contains( "_unplaced", 0 ) ) {
    int pos = chr_name.After( "_unplaced" ).Int( );

    if ( pos < 138 )
      return ( pos + 24 );
    else
      ForceAssert( 1 == 0 );
  }

  ForceAssert( 1 == 0 );
  return -1;
}



#endif
