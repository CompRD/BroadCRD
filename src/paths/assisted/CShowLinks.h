///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__ASSISTED__C_SHOW_LINKS__H
#define PATHS__ASSISTED__C_SHOW_LINKS__H

#include "Basevector.h"
#include "PairsManager.h"
#include "lookup/LookAlign.h"
#include "paths/ReadLoc.h"

/**
 * struct SLibInfo (core pairs info, as lib names, sep, devs).
 */
struct SLibInfo {
  
  SLibInfo( ) { }

  void Set( const String &pairs_file );
  
  longlong nreads_;
  vec<String> names_;
  vec<int> seps_;
  vec<int> devs_;

};

/**
 * struct SLinkInfo (info from one link, as read class, lib id).
 */
struct SLinkInfo {

  SLinkInfo( ) { }
  
  SLinkInfo( int rclass, int lib, bool rc1, bool rc2, int sep, int dev );
  
  // Generate a PrintTabular-ready "line".
  vec<String> ToLine( ) const;
  
  friend bool operator< ( const SLinkInfo &left, const SLinkInfo &right );

  int class_;
  int lib_;
  bool rc1_;
  bool rc2_;
  int sep_;
  int dev_;
  
};

/**
 * CShowLinks
 *
 * Class to interactively display linking information between two
 * oriented contigs.
 */
class CShowLinks {
  
public:
  
  CShowLinks( const vecbvec &bases,
	      const SLibInfo &ijumps,
	      const SLibInfo &iJumps,
	      const vec<look_align> &hits,
	      read_locs_on_disk &locs_file );

  // Interactive entry point.
  void GoInteractive( ) const;
  
  
private:

  // Set maps etc.
  void Setup( );

  // Find all links between signed contigs (answer is sorted).
  void FindLinks( int c1, int c2, vec<SLinkInfo> &slinks ) const;

  // Show links between the two contigs (signed, to capture orientation).
  void ShowLinks( int c1, int c2, ostream &out ) const;
  
  
private:  

  // Constant intput core data.
  const vecbvec &bases_;
  const SLibInfo &ijumps_;
  const SLibInfo &iJumps_;
  const vec<look_align> &hits_;     // may be empty
  read_locs_on_disk &locs_file_;
  
  // Maps and helpers.
  vec<int> to_hit_;                 // map contig_id to hit_id
  
};

#endif
