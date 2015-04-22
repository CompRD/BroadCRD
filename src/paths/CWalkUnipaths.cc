///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "String.h"
#include "Vec.h"
#include "VecUtilities.h"
#include "graph/Digraph.h"
#include "math/Functions.h"
#include "paths/CWalkUnipaths.h"
#include "paths/Unipath.h"

/**
 * class CWalkUnipaths
 * Constructor
 */
CWalkUnipaths::CWalkUnipaths ( const int K,
			       const vec<NormalDistribution> *cov,
			       const vec<int> *to_rc,
			       const KmerBaseBroker *kbb,
			       const vecKmerPath *unipaths,
			       const digraph *unigraph,
			       ostream *log ) :
  K_ ( K )
{
  this->SetPointers( cov, to_rc, kbb, unipaths, unigraph, log );
  this->SetDefaults( );
}

/**
 * class CWalkUnipaths
 * SetPointers
 */
void CWalkUnipaths::SetPointers( const vec<NormalDistribution> *cov,
				 const vec<int> *to_rc,
				 const KmerBaseBroker *kbb,
				 const vecKmerPath *unipaths,
				 const digraph *unigraph,
				 ostream *log )
{
  cov_ = cov;
  to_rc_ = to_rc;
  kbb_ = kbb;
  unipaths_ = unipaths;
  unigraph_ = unigraph;
  log_ = log;
}

/**
 * class CWalkUnipaths
 * SetDefaults
 */
void CWalkUnipaths::SetDefaults( )
{
  min_klen_ = 250;
  min_extensions_ = MkVec( 40*K_, 20*K_, 10*K_, K_ );
}

/**
 * class CWalkUnipaths
 * SetMinKlen
 */
void CWalkUnipaths::SetMinKlen( const int min_klen )
{
  min_klen_ = min_klen;
}

/**
 * class CWalkUnipaths
 * SetMinExtensions
 */
void CWalkUnipaths::SetMinExtensions( const vec<int> &min_extensions )
{
  min_extensions_ = min_extensions;
}

/**
 * class CWalkUnipaths
 * BuildWalks
 */
void CWalkUnipaths::BuildWalks( )
{
  this->Cleanup( );
  
  // Part 1: build walks by looking at coverages of unipaths. This
  //  step generates pairs of walks (fw and rc copies): walks_[2*ii]
  //  and walks_[1+(2*ii)] are guaranteed to be rc to each other.

  // The adjacencies.
  int n_unipaths = unipaths_->size( );
  const vec< vec<int> > &from = unigraph_->From( );
  const vec< vec<int> > &to = unigraph_->To( );

  // Sort unipaths (higher covearge on top), cli = coverage, length, id.
  if ( log_ )
    *log_ << Date( ) << ": sorting " << n_unipaths << " unipaths" << endl;
  vec< triple<float,int,int> > cli;
  cli.reserve( n_unipaths );
  for (int id=0; id<n_unipaths; id++) {
    float this_cov = (*cov_)[id].mu_;
    int this_len = (*unipaths_)[id].TotalLength( );
    cli.push_back( make_triple( this_cov, this_len, id ) );
  }
  sort( cli.rbegin( ), cli.rend( ) );

  // Keep track of placed unipaths.
  vec<int> placed( n_unipaths, -1 );   // map unipath_id -> walk_id
  walks_.reserve( n_unipaths );

  // Keep walking until no more extensions can be found.
  int last_vinit = 0;
  if ( log_ ) *log_ << Date( ) << ": building walks... " << flush;
  while( 1 ) {

    // Find first unplaced unipath.
    int vinit = -1;
    for (int ii=last_vinit; ii<cli.isize( ); ii++) {
      int uid = cli[ii].third;
      if ( placed[uid] < 0 ) {
	last_vinit = ii;
	vinit = uid;
	break;
      }
    }
    if ( vinit < 0 ) break;   // all unipaths accounted for.

    // Two walks on each iteration: fw (even walks) and rc (odd).
    IntVec walkfw;
    IntVec walkrc;

    // Walk backward.
    int uid = vinit;
    while ( uid < (int)unipaths_->size( ) ) {
      walkfw.push_back( uid );
      placed[uid] = (int)walks_.size( );

      walkrc.push_back( (*to_rc_)[uid] );   // palyndromes on fw walk only.
      if ( (*to_rc_)[uid] != uid )
	placed[ (*to_rc_)[uid]] = 1 + (int)walks_.size( );
      
      const vec<int> &prev_uids = to[uid];
      int best_id = -1;
      float best_cov = 0.0;
      for (int ii=0; ii<prev_uids.isize( ); ii++) {
	int prev_uid = prev_uids[ii];
	if ( placed[prev_uid] > -1 ) continue;
	if ( placed[ (*to_rc_)[prev_uid] ] == (int)walks_.size( ) ) continue;
	if ( best_id < 0 || (*cov_)[prev_uid].mu_ > best_cov ) {
	  best_id = prev_uid;
	  best_cov = (*cov_)[prev_uid].mu_;
	}
      }
      if ( best_id < 0 ) break;   // no next best found.
      uid = best_id;
    }
    
    // Reverse walks.
    reverse( walkfw.begin( ), walkfw.end( ) );
    reverse( walkrc.begin( ), walkrc.end( ) );

    // Now walk forward.
    uid = vinit;
    while ( uid < (int)unipaths_->size( ) ) {
      if ( uid != vinit ) {
	walkfw.push_back( uid );
	placed[uid] = (int)walks_.size( );
	
	walkrc.push_back( (*to_rc_)[uid] );   // palyndromes on fw walk only.
	if ( (*to_rc_)[uid] != uid )
	  placed[ (*to_rc_)[uid]] = 1 + (int)walks_.size( );	
      }
      
      const vec<int> &next_uids = from[uid];
      int best_id = -1;
      float best_cov = 0.0;
      for (int ii=0; ii<next_uids.isize( ); ii++) {
	int next_uid = next_uids[ii];
	if ( placed[next_uid] > -1 ) continue;
	if ( placed[ (*to_rc_)[next_uid] ] == (int)walks_.size( ) ) continue;
	if ( best_id < 0 || (*cov_)[next_uid].mu_ > best_cov ) {
	  best_id = next_uid;
	  best_cov = (*cov_)[next_uid].mu_;
	}
      }
      if ( best_id < 0 ) break;   // no next best found.
      uid = best_id;
    }
    
    // Reverse walkrc.
    reverse( walkrc.begin( ), walkrc.end( ) );

    // Add walks to list.
    walks_.push_back( walkfw );
    walks_.push_back( walkrc );

  }
  if ( log_ ) *log_ << walks_.size( ) << " walks found" << endl;
  
  // All unipaths must be accounted for.
  int n_seen = 0;
  for (int ii=0; ii<placed.isize( ); ii++)
    if ( placed[ii] > -1 ) n_seen++;
  ForceAssertEq( n_seen, n_unipaths );
  
  // Part 2: find where the walks intersect each other, and swap
  //  branches (if needed) so to generate longer walks. This addresses
  //  some sequencing-related pathologies, as leftover illumina
  //  primers that have extremely high coverage.

  // Iteratively extend walks.
  if ( log_ ) *log_ << Date( ) << ": starting extending iterations:" << endl;
  for (int pass=0; pass<min_extensions_.isize( ); pass++) {
    int n_iter = 0;
    while ( 1 ) {
      if ( log_ )
	*log_ << " iter_" << min_extensions_[pass]
	      << "." << n_iter++ << " ... " << flush;
      int n_merges = this->ExtendWalks( min_extensions_[pass] );
      if ( log_ ) *log_ << n_merges << " merges" << endl;
      if ( n_merges < 1 ) break;
    }
  }

  // Remove mate walks_.
  if ( log_ ) *log_ << Date( ) << ": removing rc walks... " << flush;
  for (int ii=0; ii<(int)walks_.size( ) / 2; ii++)
    walks_[ 1 + 2 * ii ].clear( );
  
  // Report remaining walks_.
  int n_in = walks_.size( );
  int n_out = 0;
  for (int ii=0; ii<(int)walks_.size( ); ii++)
    if ( walks_[ii].size( ) > 0 ) n_out++;
  if ( log_ ) *log_ << n_out << " walks left" << endl;

}

/**
 * class CWalkUnipaths
 * SaveContigs
 *
 * Output files:
 *   <head_out>.fastb                      contigs
 *   <head_out>.walks                      unipaths in each walk
 *   <head_out><min_klen_>.vlabels         vlabels file for the digraph below
 *   <head_out><min_klen_>.digraph.k<K_>   digraph of contigs >= min_klen_
 *   <head_out><min_klen_>.hyper.k<K_>     same, but as an HyperKmerPath
 *   <head_out><min_klen_>.fastb           contigs from the hyper above
 */
void CWalkUnipaths::SaveContigs( const String &head_out ) const
{
  // File names.
  String strK = ToString( K_ );
  String str_min_klen = ToString( min_klen_ );
  String bases_file = head_out + ".fastb";
  String walks_file = head_out + ".walks";
  String labels_file = head_out + str_min_klen + ".vlabels";
  String shaved_file = head_out + str_min_klen + ".digraph.k" + strK;
  String hyper_file = head_out + str_min_klen + ".HyperKmerPath.k" + strK;
  String hyperb_file = head_out + str_min_klen + ".fastb";
  
  // Generate contigs.
  vecbvec contigs;
  this->WalksToContigs( contigs );

  if ( log_ ) *log_ << Date( ) << ": saving contigs" << endl;
  contigs.WriteAll( bases_file );

  // Find cross points.
  if ( log_ ) *log_ << Date( ) << ": tagging cross-points" << endl;
  vec< vec<bool> > is_cross;
  this->TagCross( is_cross );
  
  // Walks.
  if ( log_ ) *log_ << Date( ) << ": saving walks" << endl;
  ofstream wout ( walks_file.c_str( ) );
  for (int contig_id=0; contig_id<(int)contigs.size( ); contig_id++) {
    const bvec &bases = contigs[contig_id];
    const IntVec &walk = walks_[contig_id];
    
    int tot_len = 0;
    for (int ii=0; ii<(int)walk.size( ); ii++)
      tot_len += (*unipaths_)[ walk[ii] ].TotalLength( );
    if ( tot_len > 0 ) tot_len += K_ - 1;
    ForceAssertEq( tot_len, (int)bases.size( ) );

    wout << "contig_" << contig_id << " (" << bases.size( ) << " bases)\n";
    if ( walk.empty( ) ) {
      wout << "\n";
      continue;
    }

    vec< vec<String> > table;
    vec<String> line = MkVec( String( "  contig" ),
			      String( "unipath_id" ),
			      String( "beg_bp" ),
			      String( "end_bp" ),
			      String( "nkmers" ),
			      String( "cov" ),
			      String( "" ) );
    table.push_back( line );
    
    int base_pos = 0;
    for (int ii=0; ii<(int)walk.size( ); ii++) {
      int klen = (*unipaths_)[ walk[ii] ].TotalLength( );
      String str_cid = "  c" + ToString( contig_id );
      String str_ii = ToString( ii ) + "/" + ToString( walk.size( ) );
      line = MkVec( str_cid + "." + str_ii,
		    ToString( walk[ii] ),
		    ToString( base_pos ),
		    ToString( base_pos + klen + K_ - 1 ),
		    ToString( klen ),
		    ToString( (*cov_)[ walk[ii] ].mu_, 1 ),
		    String( is_cross[contig_id][ii] ? "x" : "" ) );
      table.push_back( line );
      base_pos += klen;
    }
    PrintTabular( wout, table, 3, "rrrrrrl" );
    wout << "\n";
  }
  wout.close( );

  // Shaved digraph and companion labels.
  String str_min = "(contigs >= " + ToString( min_klen_ ) + " kmers)";
  if ( log_ ) *log_ << Date( ) << ": saving digraph " << str_min << endl;
  vec<String> vlabels;
  digraph shaved;
  this->ShaveGraph( shaved, vlabels );
  WRITE( labels_file, vlabels );
  BinaryWriter::writeFile( shaved_file, shaved );
  
  // HyperKmerPath of shaved graph.
  if ( log_ ) *log_ << Date( ) << ": saving HyperKmerPath " << str_min << endl;
  HyperKmerPath hyper;
  this->ShavedToHyper( shaved, hyper );
  BinaryWriter::writeFile( hyper_file, hyper );
  hyper.DumpFastb( hyperb_file, *kbb_ );

}

/**
 * class CWalkUnipaths
 * ExtendWalks
 * private
 *
 * For each walk W it tries to determine if a crossing walk W' can be
 * used to extend W on the right. Arrows are adjacencies, vertices are
 * unipaths:
 *
 *     W   ... * -> * -> * -> * -> *
 *                       |
 *                       V
 *     W'           * -> * -> * -> * -> * ...
 *
 * Only accept extensions longer than min_extension.
 */
int CWalkUnipaths::ExtendWalks( int min_extension )
{
  // Adjacencies.
  int n_unipaths = unipaths_->size( );
  const vec< vec<int> > &from = unigraph_->From( );
  const vec< vec<int> > &to = unigraph_->To( );
  
  // Map unipath_id to contig_id, and unipath_id to position_in_walk.
  vec<int> to_contig( n_unipaths, -1 );
  vec<int> to_walkpos( n_unipaths, -1 );
  for (int contig_id=0; contig_id<(int)walks_.size( ); contig_id++) {
    for (int jj=0; jj<(int)walks_[contig_id].size( ); jj++) {
      to_contig[ walks_[contig_id][jj] ] = contig_id;
      to_walkpos[ walks_[contig_id][jj] ] = jj;
    }
  }
  
  // Left/right weight (total length) for each unipath in walks_.
  vec< vec< pair<int,int> > > weights( walks_.size( ) );
  for (int contig_id=0; contig_id<(int)walks_.size( ); contig_id++) {
    const IntVec &walk = walks_[contig_id];
    vec< pair<int,int> > &weight = weights[contig_id];
    weight.reserve( walk.size( ) );
    int wleft = 0;
    int wright = 0;
    for (int ii=1; ii<(int)walk.size( ); ii++)
      wright += (*unipaths_)[ walk[ii] ].TotalLength( );
    for (int ii=0; ii<(int)walk.size( ); ii++) {
      weight.push_back( make_pair( wleft, wright ) );
      if ( ii < (int)walk.size( )-1 ) {
	wleft += (*unipaths_)[ walk[ii] ].TotalLength( );
	wright += - (*unipaths_)[ walk[ii+1] ].TotalLength( );;
      }
    }
  }
  
  // Sort contigs by length.
  vec<int> clens;
  vec< pair<int,int> > len2contig;
  clens.reserve( walks_.size( ) );
  len2contig.reserve( walks_.size( ) );
  for (int ii=0; ii<(int)walks_.size( ); ii++) {
    if ( walks_[ii].size( ) < 1 ) {
      clens.push_back( 0 );
      len2contig.push_back( make_pair( 0, ii ) );
      continue;
    }
    int clen = weights[ii][0].second+(*unipaths_)[walks_[ii][0]].TotalLength();
    clens.push_back( clen );
    len2contig.push_back( make_pair( clen, ii ) );
  }
  sort( len2contig.rbegin( ), len2contig.rend( ) );

  // Some counters.
  vec<bool> modified( walks_.size( ), false );    // modified contigs tracker
  VecIntVec extra_walks;                          // extra walks
  int n_merges = 0;                               // number of merges
  
  // Loop over all contigs (start from larger ones).
  for (int mapii=0; mapii<(int)walks_.size( ); mapii++) {
    const int contig_id = len2contig[mapii].second;
    const IntVec &walk = walks_[contig_id];
    if ( modified[contig_id] ) continue;

    // Find best candidate for extension.
    vec< triple<int,int,int> > candidates;   // wright, walkpos, id_in_from
    for (int walkpos=0; walkpos<(int)walk.size( ); walkpos++) {
      const vec<int> &froms = from[ walk[walkpos] ];
      int wleft = weights[contig_id][walkpos].first;
      int wright = weights[contig_id][walkpos].second;
      if ( wleft < 1 && wright < 1 ) continue;
      
      for (int ii=0; ii<froms.isize( ); ii++) {
	int uid = froms[ii];
	int cid = to_contig[uid];
	int wpos = to_walkpos[uid];
	
	// Already modified.
	if ( cid < 0 ) continue;
	if ( modified[cid] ) continue;

	// Crossing itself or its rc copy.
	if ( cid == contig_id ) continue;
	if ( cid == contig_id + 1 || contig_id == cid + 1 ) continue;

	// Not an usable extension (crossing walk must extend on right).
	int wR = weights[cid][wpos].second;
	if ( wR <= min_extension + wright ) continue;
	
	// Not a proper extension.
	int wL = weights[cid][wpos].first;
	if ( wL >= wleft ) continue;
      
	// A candidate
	candidates.push_back( triple<int,int,int>( wR, walkpos, ii ) );
      }
    }

    // No candidates found, leave contig.
    if ( candidates.size( ) < 1 ) continue;

    // Sort candidates (highest extension at top).
    sort( candidates.rbegin( ), candidates.rend( ) );
    
    // Perform cross over.
    n_merges++;

    const triple<int,int,int> &winner = candidates[0];
    int walkpos = winner.second;
    int uid = from[ walk[walkpos] ][winner.third];
    int cid = to_contig[uid];
    int wpos = to_walkpos[uid];
    this->CrossOver( contig_id, walkpos, cid, wpos, extra_walks, modified );
    
  } // next contig.
  
  // Add extra bits (if needed).
  if ( extra_walks.size( ) > 0 )
    copy( extra_walks.begin( ), extra_walks.end( ), back_inserter( walks_ ) );
  
  // Done.
  return n_merges;

}

/**
 * class CWalkUnipaths
 * CrossOver
 * private
 *
 * Extend first walk (contig_id) by jumping on the second walk (cid)
 * at the suture point walkpos-wpos:
 *
 *     contig_id   ... * -> * -> pos1 -> * -> *
 *                                |
 *                                V
 *     cid                  * -> pos2 -> * -> * -> * ...
 *
 */
void CWalkUnipaths::CrossOver( int contig_id,
			       int walkpos,
			       int cid,
			       int wpos,
			       VecIntVec &extra_walks,
			       vec<bool> &modified )
{
  const IntVec &walk = walks_[contig_id];
  
  // Swap chunks, deal with leftovers.
  IntVec left1 = walk;
  left1.resize( walkpos + 1 );
  
  IntVec right1 = walk;
  right1.erase( right1.begin( ), right1.begin( ) + walkpos + 1 );
  
  IntVec left2 = walks_[cid];
  left2.resize( wpos );
  
  IntVec right2 = walks_[cid];
  right2.erase( right2.begin( ), right2.begin( ) + wpos );
  
  IntVec newLong = left1;
  copy( right2.begin( ), right2.end( ), back_inserter( newLong ) );
  
  IntVec newLongRc;
  this->RcWalk( newLong, newLongRc );
  
  IntVec newShort;
  IntVec newExtra;
  if ( left2.size( ) < 1 ) newShort = right1;
  else if ( right1.size( ) < 1 ) newShort = left2;
  else {
    newShort = right1;
    newExtra = left2;
  }
  
  IntVec newShortRc;
  this->RcWalk( newShort, newShortRc );
  
  IntVec newExtraRc;
  if ( newExtra.size( ) > 0 )
    this->RcWalk( newExtra, newExtraRc );
  
  // Update walks_, and add extra bit (if needed).
  int mate_contig_id = contig_id + ( ( contig_id % 2 == 0 ) ? 1 : -1 );
  int mate_cid = cid + ( ( cid % 2 == 0 ) ? 1 : -1 );
  walks_[contig_id] = newLong;
  walks_[mate_contig_id] = newLongRc;
  walks_[cid] = newShort;
  walks_[mate_cid] = newShortRc;
  if ( newExtra.size( ) > 0 ) {
    extra_walks.push_back( newExtra );
    extra_walks.push_back( newExtraRc );
  }
  
  // Tag contigs as modified.
  modified[contig_id] = true;
  modified[mate_contig_id] = true;
  modified[cid] = true;
  modified[mate_cid] = true;
  
}

/**
 * class CWalkUnipaths
 * WalksToContigs
 * private
 */
void CWalkUnipaths::WalksToContigs( vecbvec &contigs ) const
{
  contigs.clear( );
  
  if ( log_ ) *log_ << Date( ) << ": generating contigs" << endl;
  contigs.reserve( walks_.size( ) );
  bvec empty;
  for (int ii=0; ii<(int)walks_.size( ); ii++) {
    if ( walks_[ii].empty( ) ) {
      contigs.push_back( empty );
      continue;
    }
    KmerPath pathwalk;
    for (int jj=0; jj<(int)walks_[ii].size( ); jj++)
      pathwalk.Append( (*unipaths_)[ walks_[ii][jj] ] );
    bvec cgbases = kbb_->Seq( pathwalk );
    contigs.push_back( cgbases );
  }
}

/**
 * class CWalkUnipaths
 * ContigKlen
 * private
 */
 int CWalkUnipaths::ContigKlen( int contig_id ) const
{
  int contig_len = 0;
  for (int ii=0; ii<(int)walks_[contig_id].size( ); ii++) 
    contig_len += (*unipaths_)[ walks_[contig_id][ii] ].TotalLength( );
  return contig_len;
}

/**
 * class CWalkUnipaths
 * RcWalk
 * private
 */
 void CWalkUnipaths::RcWalk( const IntVec &walk, IntVec &walk_rc ) const
{
  walk_rc.clear( );
  walk_rc.reserve( walk.size( ) );
  for (int ii=walk.size( )-1; ii>=0; ii--)
    walk_rc.push_back( (*to_rc_)[ walk[ii] ] );
}

/**
 * class CWalkUnipaths
 * TagCross
 * private
 *
 * Tag all cross points (where a "long" walk could cross-over to
 * another "long" walk, or to a different section of itself).
 */
void CWalkUnipaths::TagCross( vec< vec<bool> > &is_cross ) const
{
  int n_vertices = unigraph_->N( );

  // Default:  no cross-points.
  is_cross.clear( );
  is_cross.resize( walks_.size( ) );
  for (int ii=0; ii<(int)walks_.size( ); ii++)
    is_cross[ii].resize( walks_[ii].size( ), false );

  // Maps (unipath_id -> contig_id) and (unipath_id -> pos_in_contig).
  vec<int> to_contig( unipaths_->size( ), -1 );
  vec<int> to_cgpos( unipaths_->size( ), -1 );
  for (int contig_id=0; contig_id<(int)walks_.size( ); contig_id++) {
    if ( this->ContigKlen( contig_id ) < min_klen_ ) continue;
    for (int pos=0; pos<(int)walks_[contig_id].size( ); pos++) {
      int uid = walks_[contig_id][pos];
      to_contig[uid] = contig_id;
      to_cgpos[uid] = pos;
    }
  }
  
  // Loop over all contigs.
  for (int contig_id=0; contig_id<(int)walks_.size( ); contig_id++) {
    if ( this->ContigKlen( contig_id ) < min_klen_ ) continue;

    // Loop over all vertices in contigs' walk;
    for (int pos=0; pos<(int)walks_[contig_id].size( ); pos++) {
      int uid = walks_[contig_id][pos];

      // uid -> vid.
      const vec<int> &from = unigraph_->From( uid );
      for (int ii=0; ii<from.isize( ); ii++) {
	int vid = from[ii];
	if ( to_contig[vid] < 0 || to_cgpos[vid] < 0 ) continue;
	if ( to_contig[vid] != contig_id || to_cgpos[vid] != to_cgpos[uid]+1 ) {
	  is_cross[contig_id][pos] = true;
	  break;
	}
      }
      if ( is_cross[contig_id][pos] ) continue;
      
      // vid -> uid.
      const vec<int> &to = unigraph_->To( uid );
      for (int ii=0; ii<to.isize( ); ii++) {
	int vid = to[ii];
	if ( to_contig[vid] < 0 || to_cgpos[vid] < 0 ) continue;
	if ( to_contig[vid] != contig_id || to_cgpos[vid] != to_cgpos[uid]-1 ) {
	  is_cross[contig_id][pos] = true;
	  break;
	}
      }

    } // loop over vertices in contigs' walk.

  } // loop over all contigs.

}

/**
 * class CWalkUnipaths
 * ShaveGraph
 * private
 *
 * Shave unigraph_ (an adjacency is kept if and only if both its two
 * vertices belongs to a walk with kmer-length >= min_klen_).
 */
void CWalkUnipaths::ShaveGraph( digraph &shaved, vec<String> &vlabels ) const
{
  int n_vertices = unigraph_->N( );

  // No need to clear shaved, it will be Initalized later on.
  vlabels.clear( );
  vlabels.resize( n_vertices, "" );

  // Tag vertices to keep (and fill vlabels accordingly).
  vec<bool> to_keep( n_vertices, false );
  for (int contig_id=0; contig_id<(int)walks_.size( ); contig_id++) {
    if ( this->ContigKlen( contig_id ) < min_klen_ ) continue;
    const IntVec &walk = walks_[contig_id];
    for (int ii=0; ii<(int)walk.size( ); ii++) {
      to_keep[ walk[ii] ] = True;      
      int klen = (*unipaths_)[ walk[ii] ].TotalLength( );
      vlabels[ walk[ii] ]
	= "\\nc" + ToString( contig_id )
	+ "\\n" + ToString( klen )
	+ " kmer" + String( ( klen > 1 ? "s" : "" ) )
	+ "\\ncov: " + ToString( (*cov_)[ walk[ii] ].mu_, 1 );
    }
  }
  
  // Cut all adjacencies from/to vertices that are not tagged as keepers.
  const vec< vec<int> > &orig_froms = unigraph_->From( );
  const vec< vec<int> > &orig_tos = unigraph_->To( );
  vec< vec<int> > froms( n_vertices );
  vec< vec<int> > tos( n_vertices );
  for (int ii=0; ii<n_vertices; ii++) {
    if ( ! to_keep[ii] ) continue;
    for (int jj=0; jj<orig_froms[ii].isize( ); jj++) {
      if ( to_keep[ orig_froms[ii][jj] ] ) {
	froms[ii].push_back( orig_froms[ii][jj] );
	
      }
    }
    for (int jj=0; jj<orig_tos[ii].isize( ); jj++)
      if ( to_keep[ orig_tos[ii][jj] ] )
	tos[ii].push_back( orig_tos[ii][jj] );
  }

  // The graph.
  shaved.Initialize( froms, tos );
  
}

/**
 * class CWalkUnipaths
 * ShavedToHyper
 * private
 *
 * Convert the given digraph into an HyperKmerPath.
 */
void CWalkUnipaths::ShavedToHyper( const digraph &shaved,
				   HyperKmerPath &hyper ) const
{
  hyper.Initialize( K_, shaved, VecOfKmerPath( *unipaths_ ) );
  hyper.RemoveSmallComponents( min_klen_ );
  hyper.RemoveUnneededVertices( );
  hyper.RemoveDeadEdgeObjects( );
}

