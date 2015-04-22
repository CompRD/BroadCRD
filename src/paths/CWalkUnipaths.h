///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__C_WALK_UNIPATHS_H
#define PATHS__C_WALK_UNIPATHS_H

#include "String.h"
#include "Intvector.h"
#include "Vec.h"
#include "math/Functions.h"
#include "paths/Unipath.h"

/**
 * class CWalkUnipaths
 *
 * Class implementaing a greedy algorithm to select linear (or
 * circular) walks from a unipath graph. The idea is to "walk" on the
 * adjacency graph preferentially moving toward the highest coverage
 * unipath.
 *
 * HEURISTICS
 *
 * min_klen_: when saving contigs, save also the subdigraph (and
 *   relative vlabels) of of unigraph_ consisting of all and only the
 *   edges between vertices that belong to contigs >= min_klen_ kmers.
 *
 * min_extensions_: after walks are generated (by walking on highest
 *   copy number unipath at junctions), try to extend walks by jumping
 *   on another (extending) walk. For example (arrows are adjacencies,
 *   vertexes are unipath):
 *
 *     walk1   ... * -> * -> pos1 -> * -> *
 *                            |
 *                            V
 *     walk2            * -> pos2 -> * -> * -> * ...
 *
 *   walk1 can be extend if we jump to walk2 at pos1/pos2. A walk is
 *   not extended if the extension is too small. The code iterates on
 *   more and more stringent (smaller) values for the accepted min
 *   extension.
 */
class CWalkUnipaths {

public:

  CWalkUnipaths( const int K,
		 const vec<NormalDistribution> *cov = 0,
		 const vec<int> *to_rc               = 0,
		 const KmerBaseBroker *kbb           = 0,
		 const vecKmerPath *unipaths         = 0,
		 const digraph *unigraph             = 0,
		 ostream *log                        = 0 );
  
  // Set methods.
  void SetPointers( const vec<NormalDistribution> *cov = 0,
		    const vec<int> *to_rc               = 0,
		    const KmerBaseBroker *kbb           = 0,
		    const vecKmerPath *unipaths         = 0,
		    const digraph *unigraph             = 0,
		    ostream *log                        = 0 );
  void SetDefaults( );
  void SetMinKlen( const int min_klen );
  void SetMinExtensions( const vec<int> &min_extensions );

  // Core build method.
  void BuildWalks( );

  // Constant accessors.
  const VecIntVec &GetWalks( ) const { return walks_; }
  
  // Generate and save contigs, etc. (see .cc for a list of saved files).
  void SaveContigs( const String &head_out ) const;
  
  
private:

  void Cleanup( ) { walks_.clear( ); }
  
  // Iteratively extend walks (it returns the number of extensions).
  int ExtendWalks( int min_extension );
  
  // Extend walk contig_id by crossing over cid at walkpos-wpos.
  void CrossOver( int contig_id, int walkpos, int cid, int wpos,
		  VecIntVec &extra_walks, vec<bool> &modified );
  
  // Turn walks into contigs (empty bvecs for empty walks, to preserve ids).
  void WalksToContigs( vecbvec &contigs ) const;
  
  // Build the rc of a given walk.
  void RcWalk( const IntVec &walk, IntVec &walk_rc ) const;
  
  // Length (in kmers) of given walk.
  int ContigKlen( int contig_id ) const;

  // Find cross-points (only between "long" contigs).
  void TagCross( vec< vec<bool> > &is_cross ) const;

  // Shave graph (keep only "long" contigs, ie >= min_klen_).
  void ShaveGraph( digraph &shaved, vec<String> &vlabels ) const;
  
  // Turn shaved graph into an HyperKmerPath.
  void ShavedToHyper( const digraph &shaved, HyperKmerPath &hyper ) const;
  
  
private:

  const int K_;                           // kmer size
  const vec<NormalDistribution> *cov_;   // coverage of unipaths
  const vec<int> *to_rc_;                 // map to rc of unipths
  const KmerBaseBroker *kbb_;             // kmer base broker
  const vecKmerPath *unipaths_;           // unipaths
  const digraph *unigraph_;               // adjacency graph

  int min_klen_;                          // used to define "long" contigs
  vec<int> min_extensions_;               // minimum sizes for extensions
  
  ostream *log_;                          // log stream
  VecIntVec walks_;                       // walks
  
};

#endif
