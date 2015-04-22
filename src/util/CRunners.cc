/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "util/CRunners.h"

/**
 * CRunners
 * Constructor
 */
CRunners::CRunners( const String &PRE,
		    const String &DATA,
		    const String &RUN ) :
  PRE_ ( PRE ),
  DATA_ ( DATA ),
  RUN_ ( RUN )
{ }

/**
 * CRunners
 * Setup
 */
void CRunners::Setup( const String &SUBDIR, const String &OUTDIR )
{
  SUBDIR_ = SUBDIR;
  OUTDIR_ = OUTDIR;

  level_ = 0;
  runners_.clear( );
}

/**
 * CRunners
 * SetInOut
 */
void CRunners::SetInOut( const String &in_type, const String &out_type )
{
  in_type_ = in_type;
  out_type_ = out_type;
}

/**
 * CRunners
 * Add
 */
void CRunners::Add( const String &module, const String &args, bool last )
{
  String sub_dir = "";
  if ( in_type_ != "" ) {
    String tail = level_ == 0 ? SUBDIR_ : SUBDIR_ + "." + ToString( level_ );
    sub_dir = in_type_ + "=" + tail + " ";
  }
  String out_dir = "";
  if ( out_type_ != "" ) {
    String tail = last ? OUTDIR_ : SUBDIR_ + "." + ToString( level_ + 1 );
    out_dir = out_type_ + "=" + tail + " ";
  }
  String pdr = this->PreDataRun( ) + " ";
  String runner = module + " " + pdr + sub_dir + out_dir + args;
  
  runners_.push_back( runner );
  if ( out_type_ != "" ) level_ += 1;
}

/**
 * CRunners
 * GetRunners
 */
vec<String> CRunners::GetRunners( ) const
{
  return runners_;
}

/**
 * CRunners
 * PreDataRun
 */
String CRunners::PreDataRun( ) const
{
  return "PRE=" + PRE_ + " DATA=" + DATA_ + " RUN=" + RUN_;
}

/**
 * CRunners
 * FullRun
 */
String CRunners::FullRun( ) const
{
  return PRE_ + "/" + DATA_ + "/" + RUN_;
}

/**
 * CRunners
 * FullSub
 */
String CRunners::FullSub( bool last ) const
{
  String sub_dir = this->FullRun( ) + "/";
  sub_dir += ( last ) ? OUTDIR_ : SUBDIR_ + ToString( level_ );
  return sub_dir;
}

