/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef C_PROJECT_PARSER_H
#define C_PROJECT_PARSER_H

#include "String.h"
#include "TokenizeString.h"
#include "Vec.h"
#include "math/Functions.h"

/**
 * class CProjectParser
 *
 * It parses a given file for assembly projects dirs.
 *
 * Empty lines in the input files are ignored. Non empty lines must be
 * either comment lines (start with the reserved character "#"), or
 * lines describing an assembly (exactly one line per assembly). Empty
 * spaces are allowed within entries.
 * 
 * An assembly line is a ":" separated string with at least one entry
 * (which is interpreted as "tag"). A valid assembly line (order
 * matters!):
 *
 *  tag:pre:data:run:subdir:outdir:base:keys
 *
 * Some fields (as pre or subdir) may be omitted, as in:
 *
 *  "elephant::projects/Loxodonta:runz/work:::on_human on_dog"
 *
 * In this case pre, subdir, outdir will be set to "", and there will
 * be no keys. Notice also the allowed empty space in base.
 */
class CProjectParser {

public:

  CProjectParser( const String *PRE = 0 ) { if ( PRE ) pre_ = *PRE; }

  String Tag( ) const { return tag_; }

  String Pre( ) const { return pre_; }

  String Data( ) const { return data_; }

  String Run( ) const { return run_; }

  String Subdir( ) const { return subdir_; }

  String Outdir( ) const { return outdir_; }

  String Base( ) const { return base_; }

  String Keys( ) const { return keys_; }
  
  // Load next assembly line (return false if none is found).
  bool Load( istream &in );

  // Print detailed info (a csv one-liner of oneline = true).
  void PrintInfo( ostream &out, bool oneline = false ) const;

  // Full DATA path name.
  String FullData( ) const;

  // Full RUN path name.
  String FullRun( ) const;
  
  // Full SUBDIR path name.
  String FullSubdir( ) const;

  // Returns "PRE=pre_ DATA=data_".
  String PackagePreData( ) const;

  // Returns "PRE=pre_ DATA=data_ RUN=run_".
  String PackagePreDataRun( ) const;

  // Returns "PRE=pre_ DATA=data_ RUN=run_ SUBDIR=subdir_".
  String PackagePreDataRunSubdir( ) const;
  
  
private :
  
  String tag_;
  String pre_;
  String data_;
  String run_;
  String subdir_;
  String outdir_;
  String base_;
  String keys_;

};

/**
 * LoadVecCProjectParsers
 * 
 * Load a vector of CProjectParsers from the given file. If PRE is given,
 * set this as the PRE for all the loaded projects.
 */
vec<CProjectParser> LoadVecCProjectParsers( const String &in_file,
					    const String PRE = "" );

#endif
