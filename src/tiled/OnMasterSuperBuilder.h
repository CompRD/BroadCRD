// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef ON_MASTER_SUPER_BUILDER_H
#define ON_MASTER_SUPER_BUILDER_H

#include "AnAssemblyClass.h"
#include "ContigPair.h"
#include "Vec.h"
#include "VecInserts.h"
#include "tiled/ValidateChimpOnHuman.h"

/*
 * class om_super_builder
 *
 * Supercontig builder for the OnMaster series.
 */
class om_super_builder {
  
public:
  
  om_super_builder( );
  
  om_super_builder( const assembly *the_assembly,
		    const vec<contig_pair> *cg_pairs,
                    const map< pair<int,int>, int> *cg_pairs_idx,
		    const vec<contig_on_human> *on_master,
                    ostream *log );

  ~om_super_builder( );
  
  void SetPointers( const assembly *the_assembly,
		    const vec<contig_pair> *cg_pairs,
                    const map< pair<int,int>, int> *cg_pairs_idx,
		    const vec<contig_on_human> *on_master,
                    ostream *log );
  
  void Generate( int on_master_pos, super &the_super ) const;

  bool NextUnplaced( int &on_master_pos, int chr ) const;

  
private:

  void FindBreakingGaps( );

  bool PickNextContig( const vec<int> &contigs, 
		       map<int,int> &failures,
		       set<int> &temp_failures,
		       int &next_pos ) const;

  bool IsSalvageable( int pos ) const;
  
  bool PlaceNextContig( vec<int> &contigs,
                        vec<int> &starts,
                        vec<int> &stdevs,
			int next_pos ) const;

  bool PlaceSalvageableContig( vec<int> &contigs_in_super,
			       int pos ) const;
  
  void RefineGaps( super &the_super,
                   vec<int> &stdevs ) const;
  
private:
  
  const assembly *the_assembly_;
  const vec<contig_pair> *cg_pairs_;
  const map< pair<int,int>, int> *cg_pairs_idx_;
  const vec<contig_on_human> *on_master_;
  ostream *log_;

  vec<Bool> *placed_;
  vec<int> break_after_contig_;
};



#endif
