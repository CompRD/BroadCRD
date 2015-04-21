///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef BAM__ASSEMBLY_TO_BAMS_DEFAULTS__H
#define BAM__ASSEMBLY_TO_BAMS_DEFAULTS__H

#include "String.h"
#include "Vec.h"

/**
 * class AssemblyToBamsDefaults
 *
 * Container for various defaults used by AssemblyToBams (for example,
 * which libraries to convert). Each line contains a key descriptor,
 * followed by a list of comma separated Strings for that key (white
 * spaces are ok). For example:
 *
 *   LIBS_DIR /seq/picard
 *   LIB C009HACXX, 1, Pond-82286
 *   LIB A07GP, 1, Solexa-72798
 *   ...
 *
 * ================ LIST OF VALID KEY/ENTRIES ================
 *
 * == 0. COMMENTS ==
 *
 *   A line starting with "#" is a comment. Lines must either be
 *   comments, or start with a valid key.
 *
 * == 1. LIBRARIES ==
 *
 *   1a. LIBS_DIR: where the libs are to be found
 *   1b. LIB: each library is specified by four entries:
 *     { flowcell, lane, library, full_name }
 *
 *   Input bam for each library is loaded from the dir
 *     <libs_dir_>/<flowcell>/<DATE>/<lane>/<library>,
 *   where <DATE> is picard defined date-based directory name (there
 *   must be exactly one <DATE> in <libs_dir_>/<flowcell>).
 *
 *   If full_name is empty, then there must be exactly one .bam file
 *   in the input directory. If there are multiple .bam files, then
 *   full_name must be assigned, to specify which bam is the one to be
 *   picked.
 *
 * == 2. ASSEMBLIES ==
 *
 *   ASSEMBLY: each assembly is specified by two entries:
 *     { head_name, tag }
 *
 *   <head_name> is the full path name to head of assembly name,
 *   <tag> is a short descriptor (tag) for the assembly.
 *
 * == 3. UNIBASES ==
 *
 *   UNIBASES: each set of unibases is specified by two entries:
 *     { unibases, tag }
 *
 *   <unibases> is the full path name of the unibases' fastb,
 *   <tag> is a short descriptor (tag) for this set of unibases.
 */
class AssemblyToBamsDefaults {

public:

  AssemblyToBamsDefaults( const String &defaults_file );

  // Libraries.
  int NLibraries( ) const { return libs_.isize( ); }
  String LibBamFile( const int lib_id ) const;
  String LibName( const int lib_id, const bool brief = false ) const;

  // Assemblies.
  int NAssemblies( ) const { return assemblies_.isize( ); }
  String AssemblyHead( int ass_id ) const { return assemblies_[ass_id]; }
  String AssemblyTag( int ass_id ) const { return assemblies_tags_[ass_id]; }

  // Unibases.
  int NUnibases( ) const { return unibases_.isize( ); }
  String UnibasesFile( int uni_id ) const { return unibases_[uni_id]; }
  String UnibasesTag( int uni_id ) const { return unibases_tags_[uni_id]; }
  
  
private:

  void CleanUp( );
  
  void LoadDefaults( const String &defaults_file );

  void ParseLibsDir  ( const vec<String> &tokens );
  void ParseLib      ( const vec<String> &tokens );
  void ParseAssembly ( const vec<String> &tokens );
  void ParseUnibases ( const vec<String> &tokens );
  
  
private:

  String libs_dir_;          // where all libraries are
  vec< vec<String> > libs_;  // the libraries

  vec<String> assemblies_;       // the assemblies (head names)
  vec<String> assemblies_tags_;  // tag names for the assemblies

  vec<String> unibases_;        // unibases (full path names)
  vec<String> unibases_tags_;   // tag names for the unibases
  
};

#endif
