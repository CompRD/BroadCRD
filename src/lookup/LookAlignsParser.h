// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef LOOK_ALIGNS_PARSER_H
#define LOOK_ALIGNS_PARSER_H

#include <map>

#include "Basevector.h"
#include "Qualvector.h"
#include "ReadLocation.h"
#include "Vec.h"
#include "VecAlignmentPlus.h"
#include "lookup/LookAlign.h"



/*
 * class look_aligns_parser
 *
 * Parse vectors of (sorted!) look_align_plus'.
 */
class look_aligns_parser {
  
public:
  
  look_aligns_parser( );

  look_aligns_parser( const vec<look_align_plus> *hits,
		      const vec_alignment_plus *aligns = 0,
		      const vecbasevector *bases = 0,
		      const vecqualvector *quals = 0 );
  
  void SetPointers( const vec<look_align_plus> *hits,
		    const vec_alignment_plus *aligns = 0,
		    const vecbasevector *bases = 0,
		    const vecqualvector *quals = 0 );
  
  void PartitionHits( vec<int> &bag_begins ) const;
  
  void AppendLocs( vec<read_location> &locs ) const;

  pair<int, int> GetBracket( const vec<int> &read_ids ) const;

  
private:

  void GenerateIdToHit( ) const;
  
  int AlignId( int hit1, int hit2 ) const;
  
  bool AlignOffset( int hit1, int hit2, int &offset ) const;
  
  bool IsGoodAlign( int align_index ) const;


private:
  
  const vec<look_align_plus> *hits_;
  const vec_alignment_plus *aligns_;  // may be null
  const vecbasevector *bases_;        // may be null
  const vecqualvector *quals_;        // may be null
  mutable map<int, int> id_to_hit_;   // read_id to hit_pos map.
  
};



#endif
