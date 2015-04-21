///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__ASSISTED__C_WALK__H
#define PATHS__ASSISTED__C_WALK__H

#include "Intvector.h"
#include "paths/KmerPath.h"

/**
 * enum termination_code
 *
 * Termination codes for the iterative cloud walker.
 */
enum termination_code {
  NOT_TERMINATED,    // walk not terminated yet
  TARGET_FOUND,      // target found
  NO_EXTENSIONS,     // no more extensions found
  TOO_LONG           // extension is too long
};

/**
 * class CWalk
 *
 * It encapsulates a single walk from a given unipath. Used in the
 * recursive algorithm CloudWalk( ), in CloudWalker.{cc,h}.
 */
class CWalk {
  
public:

  CWalk( ) : termin_ ( NOT_TERMINATED ) { }
  
  CWalk( const longlong gid,
	 const KmerPath &kpath,
	 const termination_code termin = NOT_TERMINATED );
  
  void Clear( );

  void Reserve( const size_t n_extensions );

  void SetTermin( const termination_code termin ) { termin_ = termin; }
  
  // Extend walk (or create an initial walk).
  void Extend( const longlong gid,
	       const KmerPath &ext,
	       const int *over = 0,
	       const termination_code *termin = 0 );
  
  // Const accessors.
  size_t Size( ) const { return gids_.size( ); }
  const longlong &Gid( int ii ) const { return gids_[ii]; }
  const int &Gover( int ii ) const { return govers_[ii]; }
  const KmerPath &Kpath( ) const { return kpath_; }
  const termination_code &Termin( ) const { return termin_; }
  
  // Print list of gids_ (and some info), as a one liner.
  void PrintGidsInfo( const vecKmerPath &cloud_paths, ostream &out ) const;
  
  
private:
  
  LongVec gids_;             // ids in the global cloud
  IntVec govers_;            // overlaps in kmers (of size gids_.size( )-1)
  KmerPath kpath_;           // actual kmer path of walk
  termination_code termin_;  // termination flag
  
};

#endif
