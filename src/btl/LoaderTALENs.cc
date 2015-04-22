///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "Fastavector.h"
#include "String.h"
#include "TokenizeString.h"
#include "Vec.h"
#include "btl/LoaderTALENs.h"
#include "util/RunCommand.h"
// MakeDepend: dependency FastbToFinishQualb
// MakeDepend: dependency MergeFastaFilesIntoFastb

/**
 * LoadParts
 */
String LoadParts( const String &PARTS_DIR,
		  ostream *log )
{
  const String parts_fastb_file = PARTS_DIR + "/parts.fastb";
  const String parts_qualb_file = PARTS_DIR + "/parts.qualb";
  
  // Clean up.
  if ( IsRegularFile( parts_fastb_file ) ) Remove( parts_fastb_file );
  if ( IsRegularFile( parts_qualb_file ) ) Remove( parts_qualb_file );

  // Merge all into a single fastb.
  String devnull = "/dev/null";
  String theCommand
    = "MergeFastaFilesIntoFastb IN_DIR=" + PARTS_DIR
    + " OUT_DIR=" + PARTS_DIR
    + " HEAD=parts";
  RunCommandWithLog( theCommand, devnull );

  // Fake quals (all set to finished grade).
  theCommand = "FastbToFinishQualb HEAD=parts BASE_DIR=" + PARTS_DIR;
  RunCommandWithLog( theCommand, devnull );
  
  // Log parts found.
  size_t n_parts = MastervecFileObjectCount( parts_fastb_file );
  if ( log ) *log << Date( ) << ": found " << n_parts << " valid parts" << endl;

  // Return (error).
  if ( n_parts < 1 ) {
    if ( log ) *log << "\nFATAL ERROR: cannot continue with no parts\n" << endl;
    return "";
  }

  // Return (ok).
  return parts_fastb_file;
}

/**
 * LoadProducts
 */
void LoadProducts( const String &AB1_DIR,
		   const bool R2,
		   vec<String> &prods,
		   ostream &log )
{
  // All files in indir.
  vector<String> all_files = AllFiles( AB1_DIR );
  if ( all_files.size( ) < 1 ) {
    log << "Fatal error: no ab1 files found in " << AB1_DIR << "\n" << endl;
    return;
  }
    
  // Filter out non ".ab1" files.
  for (size_t ii=0; ii<all_files.size( ); ii++) {
    const String &file = all_files[ii];
    if ( ! file.Contains( ".ab1", -1 ) ) continue;
    bool fw = file.Contains( "-F1.ab1" ) || file.Contains( "-F2.ab1" );
    bool rc = file.Contains( R2 ? "-R2.ab1" : "-R1.ab1" );
    if ( ! ( fw || rc ) ) {
      log << "WARNING! Skipping " << file << " (unsupported file name)\n";
      continue;
    }

    String prod;
    if ( fw ) {
      if ( file.Contains( "-F1.ab1" ) ) prod = file.Before( "-F1.ab1" );
      else prod = file.Before( "-F2.ab1" );
    }
    else prod = file.Before( R2 ? "-R2.ab1" : "-R1.ab1" );
    
    if ( prods.size( ) < 1 || prods.back( ) != prod )
      prods.push_back( prod );
  }
  
  // Make sure all ".ab1" files are in pairs, and satisfy naming convention.
  sort( prods.begin( ), prods.end( ) );
  prods.erase( unique( prods.begin( ), prods.end( ) ), prods.end( ) );
  
  int n_valid = prods.size( );
  for (size_t ii=0; ii<prods.size( ); ii++) {
    const String &fn = prods[ii];
    
    String f1 = fn + "-F1.ab1";
    String full_f1 = AB1_DIR + "/" + f1;
    if ( ! IsRegularFile( full_f1 ) ) {
      log << "WARNING! Skipping " << fn << " (" << f1 << " not found)\n";
      n_valid--;
      continue;
    }
    
    String f2 = fn + ( R2 ? "-R2.ab1" : "-R1.ab1" );
    String full_f2 = AB1_DIR + "/" + f2;
    if ( ! IsRegularFile( full_f2 ) ) {
      log << "WARNING! Skipping " << fn << " (" << f2 << " not found)\n";
      n_valid--;
      continue;
    }      
    
    // Third read (f3) may or may not exist.
    
  }

  // Done.
  log << Date( ) << ": found " << n_valid << " sets in AB1_DIR\n" << endl;
  
}

/**
 * LoadPrimers
 */
void LoadPrimers( const String &in_fastavec_file,
		  vecbvec &primers )
{
  primers.clear( );

  vec<fastavector> as_fasta;
  LoadFromFastaFile( in_fastavec_file, as_fasta );
  primers.reserve( as_fasta.size( ) );
  for (size_t ii=0; ii<as_fasta.size( ); ii++)
    primers.push_back( as_fasta[ii].ToBasevector( ) );
}

/**
 * LoadReferences
 */
bool LoadReferences( const String &in_references_file,
		     const vec<String> &prod_names,
		     vec<String> &references,
		     ostream *log )
{
  references.clear( );
  references.resize( prod_names.size( ), "" );

  // Nothing to do
  if ( in_references_file == "skip" ) return true;
  
  // Error: no reference file found.
  if ( !IsRegularFile( in_references_file ) ) {
    if ( log ) *log << "Fatal error: references file not found.\n"
		    << "Run with REF=skip, if you do not have one.\n"
		    << endl;
    return false;
  }

  // Parse references file.
  String line;
  vec<String> tokens;
  vec<String>::const_iterator it;
  ifstream in( in_references_file.c_str( ) );
  while ( in ) {
    getline( in, line );
    if ( ! in ) break;
    Tokenize( line, ',', tokens );
    String &name = tokens[0];
    it = find( prod_names.begin( ), prod_names.end( ), name );
    if ( it == prod_names.end( ) ) continue;
    int id = distance( prod_names.begin( ), it );
    String monos = "";
    for (size_t ii=1; ii<tokens.size( ); ii++) {
      monos += tokens[ii];
      if ( ii == tokens.size( )-1 ) monos += ".t";
      else monos += " ";
    }
    references[id] = monos;
  }

  // Make sure all prod_names have a reference.
  for (size_t ii=0; ii<references.size( ); ii++) {
    if ( references[ii] == "" ) {
      if ( log ) *log << "Fatal error: no reference found for product "
		      << prod_names[ii] << "\n"
		      << endl;
      return false;
    }
  }

  // Ok.
  return true;
}

