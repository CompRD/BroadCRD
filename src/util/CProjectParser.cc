/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "String.h"
#include "TokenizeString.h"
#include "Vec.h"
#include "math/Functions.h"
#include "util/CProjectParser.h"

/**
 * CProjectParser
 * Load
 *
 * Move along the input stream until the first valid line is found. If
 * no such line is found return false. Notice that if pre_ is already
 * assigned, then its value will not be changed.
 */
bool CProjectParser::Load( istream &in )
{
  String aline;
  while ( in ) {
    getline( in, aline );
    if ( !in ) return false;
    if ( aline.Contains( "#", 0 ) ) continue;
    if ( ! aline.Contains( ":" ) ) continue;
    
    vec<String> tokens;
    vec<char> sep( 1, ':' );
    TokenizeStrictly( aline, sep, tokens );
    int ntok = tokens.size( );
      
    tag_ = tokens[0];
    if ( pre_ == "" ) pre_ = ( ntok > 1 ) ? tokens[1] : "";
    data_ = ( ntok > 2 ) ? tokens[2] : "";
    run_ = ( ntok > 3 ) ? tokens[3] : "";
    subdir_ = ( ntok > 4 ) ? tokens[4] : "";
    outdir_ = ( ntok > 5 ) ? tokens[5] :"";
    base_ = ( ntok > 6 ) ? tokens[6] : "";
    keys_ = ( ntok > 7 ) ? tokens[7] : "";

    return true;
  }

  // Should never get here.
  return false;
}
  
/**
 * CProjectParser
 * PrintInfo
 */
void CProjectParser::PrintInfo( ostream &out, bool oneline ) const
{
  String separator = oneline ? "," : "\n";

  out << "Tag=" << tag_ << separator
      << "PRE=" << pre_ << separator
      << "DATA=" << data_ << separator
      << "RUN=" << run_ << separator
      << "SUBDIR=" << subdir_ << separator
      << "OUTDIR=" << outdir_ << separator
      << "base=" << base_ << separator
      << "keys=" << keys_
      << "\n";
}

/**
 * CProjectParser
 * FullData
 */
String CProjectParser::FullData( ) const
{
  return pre_ + "/" + data_;
}

/**
 * CProjectParser
 * FullRun
 */
String CProjectParser::FullRun( ) const
{
  return this->FullData( ) + "/" + run_;
}

/**
 * CProjectParser
 * FullSubdir
 */
String CProjectParser::FullSubdir( ) const
{
  return this->FullRun( ) + "/" + subdir_;
}

/**
 * CProjectParser
 * PackagePreData
 */
String CProjectParser::PackagePreData( ) const
{
  return "PRE=" + pre_ + " DATA=" + data_;
}

/**
 * CProjectParser
 * PackagePreDataRun
 */
String CProjectParser::PackagePreDataRun( ) const
{
  return this->PackagePreData( ) + " RUN=" + run_;
}

/**
 * CProjectParser
 * PackagePreDataRunSubdir
 */
String CProjectParser::PackagePreDataRunSubdir( ) const
{
  return this->PackagePreDataRun( ) + " SUBDIR=" + subdir_;
}

/**
 * LoadVecCProjectParsers
 */
vec<CProjectParser> LoadVecCProjectParsers( const String &in_file,
					    const String PRE )
{
  ForceAssert( IsRegularFile( in_file ) );

  vec<CProjectParser> projects;
  ifstream in( in_file.c_str( ) );
  while ( in ) {
    CProjectParser newproject( &PRE );
    if ( ! newproject.Load( in ) ) break;
    projects.push_back( newproject );
  }
  in.close( );
  
  return projects;
}

