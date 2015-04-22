// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
//

// Generates indices relative to locs.

#ifndef LOCS_HANDLER_H
#define LOCS_HANDLER_H

#include "ReadLocation.h"
#include "String.h"

class lhandler {

public:

  lhandler( int n_reads, int n_contigs );

  lhandler( int n_reads, int n_contigs, const String &locs_file );

  void LoadFromFile( const String &locs_file );
  
  // For backword compatibility, set from a vec of read_locations.
  void SetFromVec( const vec<read_location> &in_locs );

  // Size of locs.
  int Size( ) const { return locs_.size( ); }

  // First loc for given contig id.
  int FirstLoc( int contig_id ) const { return flocs_[contig_id]; }

  // How many reads in the given contig.
  int ReadsInContig( int contig_id ) const { return nreads_[contig_id]; }

  // Position of read in the given contig (return -1 on error).
  int PosInContig( int read_id, int contig_id ) const;

  // How many times read_id has been placed ( >=0 ).
  int PlacedCount( int read_id ) const { return to_loc_[read_id].size( ); }
  
  // Estimate coverage (total assembled read length, total contig length).
  void Coverage( longlong &total_rlen, longlong &total_clen ) const;

  // Returns all possible placments (as loc ids) for read_id.
  vec<int> GetAllLocationIds( int read_id ) const { return to_loc_[read_id]; }
  
  // Returns location for read_id iff read is uniquely placed (else null).
  const read_location *GetPlacement( int read_id ) const;

  // Returns location at ii.
  const read_location &operator[]( int ii ) const { return locs_[ii]; }
  
  // For backward compatibility: returns vector with all locs.
  const vec<read_location> &Locs( ) const { return locs_; }

  // Returns vector with all first_locs.
  const vec<int> &FirstLocs( ) const { return flocs_; }
  
  
private:

  int n_reads_;              // needed to resize to_loc_
  int n_contigs_;            // needed to resize flocs_
  vec<read_location> locs_;  // the locs (sorted)
  vec< vec<int> > to_loc_;   // all placements for given read_id
  vec<int> flocs_;           // first loc for given contig_id
  vec<int> nreads_;          // reads count for given contig_id

};

#endif
