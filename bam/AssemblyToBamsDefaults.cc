///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "TokenizeString.h"
#include "Vec.h"
#include "VecUtilities.h"
#include "bam/AssemblyToBamsDefaults.h"
#include "system/System.h"

/**
 * AssemblyToBamsDefaults
 * Constructor
 */
AssemblyToBamsDefaults::AssemblyToBamsDefaults( const String &defaults_file )
{
  this->LoadDefaults( defaults_file );
}

/**
 * AssemblyToBamsDefaults
 * LibBamFile
 */
String AssemblyToBamsDefaults::LibBamFile( const int lib_id ) const
{
  const vec<String> &raw = libs_[lib_id];
  String flow_dir = libs_dir_ + "/" + raw[0];
  vector<String> date_dirs = AllFiles( flow_dir );
  ForceAssert( date_dirs.size( ) == 1 );
  String datedir = date_dirs[0];
  String indir = flow_dir + "/" + datedir + "/" + raw[1] + "/" + raw[2];
  String bam_file = indir;
  if ( raw.size( ) > 3 )
    bam_file += "/" + raw[3];
  else {
    vector<String> all_files = AllFiles( indir );
    vec<String> all_bams;
    for (int ii=0; ii<(int)all_files.size( ); ii++) {
      String filename = indir + "/" + all_files[ii];
      if ( filename.Contains( ".bam", -1 ) )
	all_bams.push_back( filename );
    }
    ForceAssert( all_bams.size( ) == 1 );
    ForceAssert( IsRegularFile( all_bams[0] ) );
    bam_file = all_bams[0];
  }
  return bam_file;
}

/**
 * AssemblyToBamsDefaults
 * LibName
 */
String AssemblyToBamsDefaults::LibName( const int lib_id,
					const bool brief ) const
{
  const vec<String> &lib = libs_[lib_id];
  String name = lib[0] + "." + lib[1];
  if ( ! brief ) name += "." + lib[2];
  return name;
}

/**
 * AssemblyToBamsDefaults
 * CleanUp
 * private
 */
void AssemblyToBamsDefaults::CleanUp( )
{
  libs_.clear( );

  assemblies_.clear( );
  assemblies_tags_.clear( );

  unibases_.clear( );
  unibases_tags_.clear( );
}

/**
 * AssemblyToBamsDefaults
 * LoadDefaults
 * private
 */
void AssemblyToBamsDefaults::LoadDefaults( const String &defaults_file )
{
  this->CleanUp( );
  
  String line;
  vec<String> tokens;
  vec<char> seps = MkVec( ' ', '\t', ',' );
  
  ifstream in( defaults_file.c_str( ) );
  while ( in ) {
    getline( in, line );
    if ( ! in ) break;
    if ( line.Contains( "#", 0 ) ) continue;
    Tokenize( line, seps, tokens );
    if ( tokens.size( ) < 1 ) continue;
    if ( tokens[0] == "LIBS_DIR" ) {
      this->ParseLibsDir( tokens );
      continue;
    }
    if ( tokens[0] == "LIB" ) {
      this->ParseLib( tokens );
      continue;
    }
    if ( tokens[0] == "ASSEMBLY" ) {
      this->ParseAssembly( tokens );
      continue;
    }
    if ( tokens[0] == "UNIBASES" ) {
      this->ParseUnibases( tokens );
      continue;
    }
    ForceAssert( 1 == 0 );   // not a valid line
  }
  in.close( );
}

/**
 * AssemblyToBamsDefaults
 * ParseLibsDir
 * private
 */
void AssemblyToBamsDefaults::ParseLibsDir( const vec<String> &tokens )
{
  ForceAssert( tokens.size( ) > 1 );
  libs_dir_ = tokens[1];
}

/**
 * AssemblyToBamsDefaults
 * ParseLib
 * private
 */
void AssemblyToBamsDefaults::ParseLib( const vec<String> &tokens )
{
  ForceAssert( tokens.size( ) > 3 );
  vec<String> lib = MkVec( tokens[1], tokens[2], tokens[3] );
  if ( tokens.size( ) > 4 ) lib.push_back( tokens[4] );
  libs_.push_back( lib );
}

/**
 * AssemblyToBamsDefaults
 * ParseAssembly
 * private
 */
void AssemblyToBamsDefaults::ParseAssembly( const vec<String> &tokens )
{
  ForceAssert( tokens.size( ) > 2 );
  assemblies_.push_back( tokens[1] );
  assemblies_tags_.push_back( tokens[2] );
}

/**
 * AssemblyToBamsDefaults
 * ParseUnibases
 * private
 */
void AssemblyToBamsDefaults::ParseUnibases( const vec<String> &tokens )
{
  ForceAssert( tokens.size( ) > 2 );
  unibases_.push_back( tokens[1] );
  unibases_tags_.push_back( tokens[2] );
}
