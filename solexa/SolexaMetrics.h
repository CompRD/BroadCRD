/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/// Define an API for storing and retrieving lane-level Solexa metrics.  Let db be a
/// solexa_metric_db.  Let m be a metric name (in double quotes).
///
/// db.Defined(m)                  - Determine if value of m is defined.
///
/// db.Value(m)                    - Return value of m, which is a String.
///                                  Assert if not defined.
///
/// db.GetValue( m, val )          - Put value of m in val or return false.
///
/// db.SetValue( m, val )          - Set m to val.
///
/// db(fn)                         - Initialize db from file fn.
///
/// db.WriteMetrics(fn, append)    - Write db to file fn, possibly appending.
///
/// Every metric in use should be defined in the human-readable file
/// "metric_defs".
///
/// Lane-level metrics are to be stored in a file FC.LANE.metrics.
///
/// Note that this is completely general - it has nothing to do with Solexa or
/// even with sequencing.
///
/// db.Value(m) will also accept pseudo-metrics of the form m[i].  If m is a metric
/// and its value has the form {x0,...,xn}, then m[i] is treated as having value xi.
/// This also works for m[i][j].

#ifndef SOLEXA_METRICS_H
#define SOLEXA_METRICS_H

#include "CoreTools.h"
#include "CommonSemanticTypes.h"
#include "ParseSet.h"
#include <map>

/// Definition of general macros.

class solexa_metric_db {

public:

  solexa_metric_db( ) { }

  /// Read metrics file from disk.  Not an error if the file does not
  /// exist; you get an empty database.  This allows updating of
  /// metrics file by many programs without having to be careful about
  /// who goes first.
  solexa_metric_db( const filename_t& f, bool requireFileAlreadyExists=true );

  bool Defined( const String& key ) const;

  String Value( const String& key ) const;

  int ValueInt( const String& key ) const
  {    return Value(key).Int( );    }

  double ValueDouble( const String& key ) const
  {    return Value(key).Double( );    }

  bool ValueBool( const String& key ) const
  {    String v = Value(key);
  if ( v == "true" ) return true;
  if ( v == "false" ) return false;
  cout << "ValueBool: illegal value \"" << v << "\" for metric "
       << key << endl;
  TracebackThisProcess( );  return false;  }

  vec<double> ValueVecDouble( const String& key ) const
  { vec<double> result; ParseDoubleSet(key,result,false); return result; }

  vec< vec<double> > ValueVecVecDouble( const String& key ) const;

  bool GetValue( const String& key, String& val ) const;

  void SetValue( const String& key, const String& val );


  /// This should work with any numerical type.
  template<class T>
  void SetValue( const String& key, const T val )
  {    SetValue( key, ToString(val) );    }

  /// Unimplemented: avoid possible confusion between int and T *.
  template<class T>
  void SetValue( const String& key, const T * val );

  /// Need this implementation to avoid confusion with templated version.
  void SetValue( const String& key, const char * val )
  {    SetValue(key, String(val)); }

  /// Need this implementation to avoid confusion with templated version.
  void SetValue( const String& key, char * val )
  {    SetValue(key, String(val)); }

  void SetValueBool( const String& key, const bool val )
  {    SetValue( key, ( val ? "true" : "false" ) );    }

  void ReadMetrics( const filename_t& f );

  void WriteMetrics( const filename_t& f, Bool append = False ) const;

  void WriteMetrics( ostream& out ) const;


  static void WriteMetricsMult( const filename_t& f, const vec< solexa_metric_db >& dbs );
  static void ReadMetricsMult( const filename_t& f, vec< solexa_metric_db >& dbs );

private:

  void ReadMetrics( ifstream& in, const filename_t& f, int maxLines = -1 );

  typedef map<String,String> maptype;
  maptype m_;

};  // class solexa_metric_db

#endif
