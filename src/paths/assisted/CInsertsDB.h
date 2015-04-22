///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__ASSISTED__C_INSERTS_DB_H
#define PATHS__ASSISTED__C_INSERTS_DB_H

#include "SeqInterval.h"
#include "SupersHandler.h"
#include "Vec.h"
#include "paths/ReadLoc.h"

/**
 * class CInsertsDB
 *
 * Deal with all the inserts in a super (providing localization for
 * inserts). Note that it only stores the fw loc from each valid pair
 * in one super.
 *
 * WARNING: the method SetSuper( ) contains a omp parallel for loop.
 */
class CInsertsDB {

public:
  
  CInsertsDB( const int MIN_SEP = 50,
	      const int MAX_SEP = 15000,
	      const int *super_id = 0,
	      const shandler *supers = 0,
	      read_locs_on_disk *locs_parser = 0 );
  
  void SetPointers( const shandler *supers, read_locs_on_disk *locs_parser );
  
  // Set super id, and fill core data structures (sis_, locs_).
  void SetSuper( const int *super_id );
  
  // Core accessors (for convenience, Size( ) clones NLocsTotal( )).
  size_t Size( ) const { return locs_.size( ); }
  const seq_interval &SeqInterval( int ii ) const { return sis_[ii]; }
  const read_loc &ReadLoc( int ii ) const { return locs_[ii]; }

  // Type of loc.
  String LocTypeAsString( int ii ) const;
  int LocTypeAsInt( int ii ) const { return sis_[ii].IntervalId( ); }
  bool IsValidType( int ii ) const { return sis_[ii].IntervalId( ) == 0; }
  bool IsSeparatedType( int ii ) const { return sis_[ii].IntervalId( ) == 1; }
  bool IsLonerType( int ii ) const { return sis_[ii].IntervalId( ) == 2; }

  // Loc counts (total, and broken by type, see .cc for details).
  int NLocsTotal( ) const;
  int NLocsOfType( int loc_type ) const;
  
  // Ids/Pos of contigs containing insert ends (see .cc for details).
  pair<int,int> ContigIds( const int ii ) const;
  pair<int,int> ContigPos( const int ii ) const;

  // Map contig_pos to insert_id (see .cc for details).
  vec< triple<int,int,int> > CposToInserts( ) const;

  // PrintTabular-ready lines with insert's info (see .cc for details).
  vec<String> LineInfo( const int ii ) const;
  vec<String> LineInfoAlt( const int ii ) const;

  // Cloud of localized inserts for given gap (see .cc for details).
  void BuildCloud( const int gap_id,
		   vec<int> &si_ids,
		   vec<int> *r_ids = 0 ) const;
  
  
private:
  
  // Decide if pair is usable, and in case return a valid seq_interval.
  bool IsValidFw( const read_loc &loc, seq_interval &si ) const;
  bool IsSeparated( const read_loc &loc, seq_interval &si ) const;
  bool IsLoner( const read_loc &loc, seq_interval &si ) const;
  
  
private:

  // Heuristics.
  int MIN_SEP_;
  int MAX_SEP_;
  
  // Const pointers to external data.
  const shandler *supers_;
  read_locs_on_disk *locs_parser_;

  // Local data (for one super).
  int super_id_;
  vec<read_loc> locs_;       // one loc per pair
  vec<seq_interval> sis_;    // interval_id encodes loc type
  
};

#endif
