/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
#include <fstream>
#include <map>
#include "String.h"
#include "system/Types.h"

#ifndef __INCLUDE_SYSTEM_PROPERTIES_H
#define __INCLUDE_SYSTEM_PROPERTIES_H

/**
   Class: Properties

   A set of (attribute,value) pairs, serializable to a file.
*/
class Properties: public map<String,String> {
 public:
  typedef map<String,String> PARENT;

  Properties() { }
  Properties( const String& fname ) { Load( fname ); }
  Properties( const Properties& p ) : PARENT(p) { }
  Properties& operator=( const Properties& p ) {
    PARENT::operator=( p );
    return *this;
  }

  void SetBool( const String& prop, Bool val ) {
    (*this)[prop] = val ? "True" : "False";
  }

  Bool GetBool( const String& prop ) {
    String val = (*this)[prop];
    ForceAssert( val == "True" || val == "False" );
    return val == "True";
  }

  void SetDouble( const String& prop, double val ) {
    char buf[128];
    sprintf( buf, "%f", val );
    (*this)[prop] = String(buf);
  }

  double GetDouble( const String& prop ) {
    String val = (*this)[prop];
    double valDbl;
    if ( sscanf( val.c_str(), "%lf", &valDbl ) != 1 )
      ForceAssert( 0 );
    return valDbl;
  }

  Bool IsDefined( const String& prop ) { return this->find( prop ) != this->end(); }
  
  void Store( const String& fname ) {
    ofstream f( fname.c_str() );
    for ( const_iterator it = begin(); it != end(); it++ )
      f << it->first << endl << (it->second == "" ? "NULL" : it->second) << endl;
  }

  void Load( const String& fname ) {
    clear();
    ifstream f( fname.c_str() );
    while ( true ) {
      String attr, val;
      f >> attr;
      if ( !f || attr == "") break;
      f >> val;
      if ( val == "NULL" ) val = "";
      (*this)[attr] = val;
    }
  }
};  // class Properties

#endif
// #ifndef __INCLUDE_SYSTEM_PROPERTIES_H
