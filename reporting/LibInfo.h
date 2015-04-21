// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
//

#ifndef LIB_INFO_H
#define LIB_INFO_H

#include "String.h"



/*
 * lib_info
 *
 * Compact container for the core data of a library (name, how many
 * reads in library, etc.)
 */
class lib_info {

public:

  lib_info( ) :
    name_ ( "" ), count_ ( 0 ), size_ ( 0 ), stdev_ ( 0 ) { }
  
  lib_info( const String &name, int size, int stdev ) :
    count_ ( 0 ) { name_ = name; size_ = size; stdev_ = stdev;
  }
  
  void Set( const String &name, int size, int stdev ) {
    count_ = 0; name_ = name; size_ = size; stdev_ = stdev;
  }
  
  void SetName( const String &name ) { name_ = name; }

  void SetCount( int count ) { count_ = count; }

  void SetSize( int size ) { size_ = size; }
  
  void SetStdev( int stdev ) { stdev_ = stdev; }

  void Increment ( longlong amount = 1 ) { count_ += amount; }
  String Name( ) const { return name_; }

  longlong Count( ) const { return count_; }

  int Size( ) const { return size_; }

  int Stdev( ) const { return stdev_; }
  
  String PrettySize( ) const {
    if ( size_ % 1000 == 0 )
      return ( ToString( size_ / 1000 ) + "Kb" );
    return ( ToString( float( size_ ) / 1000.0, 1 ) + "Kb" );
  }

  void PrettyPrint( ostream &out, bool new_line = true ) const {
    out << name_ << " ("
	<< size_ << " +/- "
	<< stdev_ << ") "
	<< count_ << ( new_line ? "\n" : "" );
  }
  
  friend bool operator== ( const lib_info &left, const lib_info &right ) {
    return ( left.name_ == right.name_ );
  }
  
  friend bool operator< ( const lib_info &left, const lib_info &right ) {
    return ( left.name_ < right.name_ );
  }
  
  
private:

  String name_;    // library name
  longlong count_; // how many reads (or bases) belong to this library
  int size_;       // insert size
  int stdev_;      // insert stdev

};



/*
 * order_Size_Kb
 *
 * This only works if name_ is of the type nKb (where n is an integer):
 * sort by library size ( the "n" in nKb).
 */
struct order_Size_Kb :
  public binary_function<const lib_info &, lib_info &, bool>
{
public:
  bool operator() ( const lib_info &left, const lib_info &right ) {
    int left_size = left.Name( ).Before( "Kb" ).Int( );
    int right_size = right.Name( ).Before( "Kb" ).Int( );
    return ( left_size < right_size );
  }
};



#endif
