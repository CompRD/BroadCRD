///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "Fastavector.h"
#include "Superb.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"
#include "paths/KmerPathInterval.h"
#include "paths/SaveScaffoldGraph.h"
#include "paths/assisted/CInsertsDB.h"
#include "paths/assisted/CloudWalker.h"
#include "paths/assisted/CWalk.h"
#include "paths/assisted/KmerPathPerfectOverlaps.h"

/**
 * OverlappingBuddies
 */
void OverlappingBuddies( const SCloud &scloud,
			 const KmerPath &kpath,
			 vec<int> &buddies )
{
  buddies.clear( );
  
  const vecKmerPath &cloud_paths = scloud.paths_;
  const vec<tagged_rpint> &cloud_pathsdb = scloud.pathsdb_;
  
  // First and last kmer in kpath.
  for (int pass=0; pass<2; pass++) {
    longlong kid = 0;
    unsigned short kpos = 0;
    if ( pass == 0 ) {
      kpos = 0;
      kid = kpath.FirstSegment( ).Start( );
    }
    else {
      kpos = kpath.NSegments( ) - 1;
      kid = kpath.LastSegment( ).Stop( );
    }
    
    // Build list of candidates locs.
    vec<longlong> locs;
    Contains( cloud_pathsdb, kid, locs );
    for (int jj=0; jj<locs.isize( ); jj++) {
      const tagged_rpint &t_rpi = cloud_pathsdb[ locs[jj] ];
      if ( t_rpi.Rc( ) ) continue;
      unsigned short rpos = t_rpi.PathPos( );
      longlong read_id = t_rpi.ReadId( );
      const KmerPath &rpath = cloud_paths[read_id];
      
      // Only accept perfect overlaps.
      if ( PerfectOverlap( kpath, rpath, kpos, rpos ) )
	buddies.push_back( t_rpi.ReadId( ) );
    }
  }
  
  // Remove duplicates.
  sort( buddies.begin( ), buddies.end( ) );
  buddies.erase( unique( buddies.begin( ), buddies.end( ) ), buddies.end( ) );
  
}

/**
 * SecondaryCloud
 */
int SecondaryCloud( const int level,
		    const SCloud &scloud,
		    vec< pair<int,int> > &ids2dist )
{
  const vecKmerPath &cloud_paths = scloud.paths_;
  const vec<tagged_rpint> &cloud_pathsdb = scloud.pathsdb_;
  
  // Sort ids.
  const size_t n_ids_init = ids2dist.size( );
  if ( ! is_sorted( ids2dist.begin( ), ids2dist.end( ) ) )
    sort( ids2dist.begin( ), ids2dist.end( ) );
  
  // Find candidates.
  vec<int> loc_ids;
  for (size_t ii=0; ii<n_ids_init; ii++) {
    const KmerPath &kpath = cloud_paths[ ids2dist[ii].first ];

    vec<int> buddies;
    OverlappingBuddies( scloud, kpath, buddies );
    copy( buddies.begin( ), buddies.end( ), back_inserter( loc_ids ) );
  }

  // Remove duplicates.
  sort( loc_ids.begin( ), loc_ids.end( ) );
  loc_ids.erase( unique( loc_ids.begin( ), loc_ids.end( ) ), loc_ids.end( ) );
  
  // Select ids to add.
  const size_t n_loc = loc_ids.size( );
  vec< pair<int,int> >::iterator it;
  vec<bool> select( n_loc, false );
  for (size_t ii=0; ii<n_loc; ii++) {
    pair<int,int> finder = make_pair( loc_ids[ii], -1 );
    it = lower_bound( ids2dist.begin( ), ids2dist.end( ), finder );
    if ( it == ids2dist.end( ) || it->first != loc_ids[ii] )
      select[ii] = true;
  }

  int n_sel = 0;
  for (size_t ii=0; ii<select.size( ); ii++)
    if ( select[ii] ) n_sel++;
  
  // Add ids.
  ids2dist.reserve( ids2dist.size( ) + size_t( n_sel ) );
  for (size_t ii=0; ii<select.size( ); ii++)
    if ( select[ii] ) ids2dist.push_back( make_pair( loc_ids[ii], level ) );
  
  // Done.
  return n_sel;
  
}

/**
 * LocalizedGlobalCloud
 */
void LocalizedGlobalCloud( const int gap_id,
			   const int gcid1,
			   const int gcid2,
			   const int MAX_DIST,
			   const SCloud &scloud,
			   const CInsertsDB &ins_db,
			   const vecKmerPath &jpaths,
			   vec< pair<int,int> > &ids2dist )
{
  ids2dist.clear( );

  const vecKmerPath &cloud_paths = scloud.paths_;
  const vec<tagged_rpint> &cloud_pathsdb = scloud.pathsdb_;

  // Localized inserts (both insert ids, and ids of end-reads).
  vec<int> local_pairs;
  vec<int> local_pairs_read_ids;
  ins_db.BuildCloud( gap_id, local_pairs, &local_pairs_read_ids );
  
  // Ids of flanking contigs are always in.
  vec<int> loc_ids = MkVec( gcid1, gcid2 );

  // Zero level ids (directly touching a jump read).
  for (int ii=0; ii<local_pairs_read_ids.isize( ); ii++) {
    const int sig_id = local_pairs_read_ids[ii];
    const int jump_id = sig_id < 0 ? - sig_id - 1 : sig_id;
    KmerPath rc_path;
    if ( sig_id < 0 ) {
      rc_path = jpaths[jump_id];
      rc_path.Reverse( );
    }
    const KmerPath &jpath = ( sig_id < 0 ) ? rc_path : jpaths[jump_id];

    vec<int> buddies;
    OverlappingBuddies( scloud, jpath, buddies );
    copy( buddies.begin( ), buddies.end( ), back_inserter( loc_ids ) );
  }

  // Remove duplicates.
  sort( loc_ids.begin( ), loc_ids.end( ) );
  loc_ids.erase( unique( loc_ids.begin( ), loc_ids.end( ) ), loc_ids.end( ) );

  // Initial setup of ids2dist (zero level sequences).
  ids2dist.reserve( loc_ids.size( ) );
  for (size_t ii=0; ii<loc_ids.size( ); ii++)
    ids2dist.push_back( make_pair( loc_ids[ii], 0 ) );

  // Add secondary cloud, up to MAX_DIST.
  for (int dist=1; dist<=MAX_DIST; dist++)
    if ( SecondaryCloud( dist, scloud, ids2dist ) < 1 )
      break;
  
}

/**
 * AllExtensions
 */
void AllExtensions( const KmerPath &read,
		    const SCloud &scloud,
		    const vec<int> &select,
		    vec<longlong> &extensions,
		    vec<KmerPath> &extended_paths )
{
  extensions.clear( );
  extended_paths.clear( );

  const vec<int> &to_rc = scloud.to_rc_;
  const vecKmerPath &cloud_paths = scloud.paths_;
  const vec<tagged_rpint> &cloud_pathsdb = scloud.pathsdb_;
  
  // Find candidate locs.
  vec<longlong> locs;
  const longlong kpi = read.LastSegment( ).Stop( );
  Contains( cloud_pathsdb, kpi, locs );
  
  // Parse locs, adding extensions when found.
  for (int jj=0; jj<locs.isize( ); jj++) {
    const tagged_rpint &t_rpi = cloud_pathsdb[locs[jj]];

    // Only fw aligns.
    if ( t_rpi.Rc( ) ) continue;
    
    // Not in the local cloud.
    int cid = t_rpi.ReadId( );
    if ( ! binary_search( select.begin( ), select.end( ), cid ) ) continue;
    
    // Not a proper extension.
    int cpos = t_rpi.PathPos( );
    int rpos = read.NSegments( ) - 1;
    const KmerPath &cpath = cloud_paths[cid];
    if ( ! IsProperExtension( cpath, read, cpos, rpos ) ) continue;
    
    // Extension found.
    KmerPath kp_ext = GlueKmerPaths( read, cpath, rpos, cpos );
    
    // This extension is a shorter version of an already found extension.
    bool subsumed = false;
    for (int kk=0; kk<extended_paths.isize( ); kk++) {
      const KmerPath &epath = extended_paths[kk];
      if ( IsReadEmbedded( kp_ext, epath, 0, 0 ) ) {
	subsumed = true;
	break;
      }
    }
    if ( subsumed ) continue;

    // This extension replaces a shorter already found extension.
    int replaces = -1;
    for (int kk=0; kk<extended_paths.isize( ); kk++) {
      const KmerPath &epath = extended_paths[kk];
      if ( IsProperExtension( kp_ext, epath, 0, 0 ) ) {
	replaces = kk;
	break;
      }
    }
    
    // Add extension.
    if ( replaces < 0 ) {
      extensions.push_back( locs[jj] );
      extended_paths.push_back( kp_ext );
    }
    else {
      extensions[replaces] = locs[jj];
      extended_paths[replaces] = kp_ext;
    }

  } // loop over all candidate locs.

}

/**
 * CloudWalk
 */
void CloudWalk( const int max_klen,
		const int max_depth,
		const int target_id,
		const SCloud &scloud,
		const vec<int> &select,		
		int depth,
		vec<CWalk> &walks,
		ostream *log )
{
  const vecKmerPath &cloud_paths = scloud.paths_;
  const vec<tagged_rpint> &cloud_pathsdb = scloud.pathsdb_;
  const int n_walks = walks.size( );

  // Max depth reached.
  depth++;
  String str_depth = "[" + ToString( depth ) + "] ";
  if ( depth >= max_depth ) {
    if ( log ) *log << str_depth << " max depth reached, stopping here\n";
    for (size_t ii=0; ii<walks.size( ); ii++) {
      if ( walks[ii].Termin( ) == NOT_TERMINATED ) {
	walks[ii].SetTermin( TOO_LONG );
	if ( log ) {
	  *log << str_depth << " " << ii << "." << n_walks << "   ";
	  walks[ii].PrintGidsInfo( cloud_paths, *log );
	}
      }
    }
    return;
  }

  if ( log ) *log << str_depth << "Starting CloudWalk\n";
   
  // This loop may cause walks to increase in size.
  for (int ii=0; ii<n_walks; ii++) {
    if ( walks[ii].Termin( ) != NOT_TERMINATED ) {
      if ( log ) {
	*log << str_depth << " " << ii << "." << n_walks << "   ";
	walks[ii].PrintGidsInfo( cloud_paths, *log );
      }
      continue;
    }
     
    // Find all extensions.
    if ( log ) *log << str_depth << "extending walk " << ii << "\n";
     
    vec<longlong> ext_ids;
    vec<KmerPath> ext_paths;
    CWalk init_walk = walks[ii];   // a copy (walks[ii] can change).
    const KmerPath &init_kpath = init_walk.Kpath( );
    const int init_len = init_kpath.TotalLength( );
    AllExtensions( init_kpath, scloud, select, ext_ids, ext_paths );
    
    // No extensions found.
    if ( ext_ids.size( ) < 1 ) {
      walks[ii].SetTermin( NO_EXTENSIONS );
      if ( log ) {
	*log << str_depth << " " << ii << "." << n_walks << "   ";
	walks[ii].PrintGidsInfo( cloud_paths, *log );
      }
      continue;
    }
    
    // Add extensions.
    for (int jj=0; jj<ext_paths.isize( ); jj++) {
      const longlong &id_ext = cloud_pathsdb[ext_ids[jj]].ReadId( );
      const KmerPath &kp_ext = ext_paths[jj];

      // Compute overlap.
      int jjcloud_len = cloud_paths[id_ext].TotalLength( );
      int ext_len = kp_ext.TotalLength( );
      int overlap = init_len + jjcloud_len - ext_len;
      ForceAssert( ext_len > init_len );
       
      // The first extension replaces walks[ii], the next are added to walks.
      if ( jj > 0 ) walks.push_back( init_walk );
      CWalk &current_walk = ( jj == 0 ) ? walks[ii] : walks[walks.size( )-1];

      // Extend walk.
      current_walk.Extend( id_ext, kp_ext, &overlap );
      if ( log ) {
	String ifnew = ( jj > 0 ? " new" : "" );
	*log << str_depth << " " << ii << "." << n_walks << ifnew << "   ";
	current_walk.PrintGidsInfo( cloud_paths, *log );
      }

      // Check termination flag.
      if ( id_ext == target_id ) {
	current_walk.SetTermin( TARGET_FOUND );
	if ( log ) *log << str_depth << "TARGET_FOUND\n";
	continue;
      }
      if ( kp_ext.KmerCount( ) > max_klen ) {
	current_walk.SetTermin( TOO_LONG );
	if ( log ) *log << str_depth << "TOO_LONG\n";
	continue;
      }
    }
    
    // Recursion.
    CloudWalk( max_klen, max_depth, target_id, scloud,
	       select, depth, walks, log );
    
  } // loop over all walks.

}

/**
 * CloudWalkWrapper
 */
void CloudWalkWrapper( const int u1,
		       const int u2,
		       const int max_klen,
		       const int max_depth,
		       const bool verbose_log,
		       const SCloud &scloud,		       
		       const vec<int> &local_ids,
		       vec<CWalk> &walks,
		       ostream *log )
{
  walks.clear( );
  
  ostream *vlog = verbose_log ? log : 0;

  const vecKmerPath &cloud_paths = scloud.paths_;
  const vec<tagged_rpint> &cloud_pathsdb = scloud.pathsdb_;
  const int u1len = cloud_paths[u1].KmerCount( );
  const int u2len = cloud_paths[u2].KmerCount( );
  
  int effective_max_klen = max_klen + u1len + u2len;

  if ( u1 == u2 ) {
    walks.push_back( CWalk( u1, cloud_paths[u1], TARGET_FOUND ) );

    if ( log )
      *log << "No need to cloud walk, insert is embedded in " << u1
	   << " (" << u2len << " kmers)\n"
	   << "\n";
    
    return;
  }
  
  if ( log )
    *log << "Cloud walking " << u1
	 << " (" << u1len << " kmers) "
	 << " to " << u2
	 << " (" << u2len << " kmers)\n";
  
  int depth = 0;   // keep track of depth (recursion level)
  walks.resize( 1, CWalk( u1, cloud_paths[u1] ) );
  CloudWalk( effective_max_klen, max_depth, u2, scloud, local_ids,
	     depth, walks, vlog );
  
  if ( log ) {
    *log << " Found " << walks.size( ) << " walks:";
    vec<int> n_terms( 4, 0 );
    for (int jj=0; jj<walks.isize( ); jj++)
      n_terms[ walks[jj].Termin( ) ]++;
    ForceAssertEq( n_terms[0], 0 );
    if ( n_terms[1] > 0 ) *log << "   " << n_terms[1] << " full";
    if ( n_terms[2] > 0 ) *log << "   " << n_terms[2] << " partial";
    if ( n_terms[3] > 0 ) *log << "   " << n_terms[3] << " too long";
    *log << "\n";
  }
  
}

/**
 * AssembleCloud
 */
void AssembleCloud( const int max_klen,
		    const int max_depth,
		    const bool verbose_log,
		    const int u1,
		    const int u2,
		    const SCloud &scloud,
		    const vec<int> &local_ids,
		    const CInsertsDB &ins_db,
		    vecKmerPath &closures,
		    vec<closure_t> &closure_types,
		    ostream *log )
{
  closures.clear( );

  const vec<int> &to_rc = scloud.to_rc_;
  
  // Try to walk from u1 to u2.
  vec<CWalk> walks;
  CloudWalkWrapper( u1, u2, max_klen, max_depth, verbose_log,
		    scloud, local_ids, walks, log );
  
  int n_closed = 0;
  for (int ii=0; ii<walks.isize( ); ii++)
    if ( walks[ii].Termin( ) == TARGET_FOUND )
      n_closed++;
  
  // Full closure(s) found, stopping here (is this really correct?).
  if ( n_closed > 0 ) {
    closures.reserve( n_closed );
    for (int ii=0; ii<walks.isize( ); ii++) {
      if ( walks[ii].Termin( ) == TARGET_FOUND ) {
	closures.push_back( walks[ii].Kpath( ) );
	closure_types.push_back( FULL_CLOSURE );
      }
    }
    return;
  }
  
  closures.reserve( walks.size( ) );
  for (int ii=0; ii<walks.isize( ); ii++) {
    closures.push_back( walks[ii].Kpath( ) );
    closure_types.push_back( LEFT_CLOSURE );
  }
  
  // Try to walk from u2 to u1 (no closures found, add partial closures).
  int u2rc = to_rc[u2];
  int u1rc = to_rc[u1];      
  CloudWalkWrapper( u2rc, u1rc, max_klen, max_depth, verbose_log,
		    scloud, local_ids, walks, log );
  
  closures.reserve( closures.size( ) + walks.size( ) );
  for (int ii=0; ii<walks.isize( ); ii++) {
    KmerPath rc_walk = walks[ii].Kpath( );
    rc_walk.Reverse( );
    closures.push_back( rc_walk );
    closure_types.push_back( RIGHT_CLOSURE );
  }
}

/**
 * AssembleCloud
 */
void AssembleCloud( const int max_klen,
		    const int max_depth,
		    const bool verbose_log,
		    const int insert_id,
		    const SCloud &scloud,
		    const vec<int> &local_ids,
		    const CInsertsDB &ins_db,
		    const vecKmerPath &jpaths,
		    vecKmerPath &closures,
		    vec<closure_t> &closure_types,
		    ostream *log )
{
  closures.clear( );
  closure_types.clear( );

  const vec<int> &to_rc = scloud.to_rc_;
  const vecKmerPath &cloud_paths = scloud.paths_;
  const vec<tagged_rpint> &cloud_pathsdb = scloud.pathsdb_;

  // Prime the walks.
  int id1 = ins_db.ReadLoc( insert_id ).ReadId( );
  vec<longlong> prim1;
  vec<KmerPath> ext1;
  AllExtensions( jpaths[id1], scloud, local_ids, prim1, ext1 );

  int id2 = ins_db.ReadLoc( insert_id ).PartnerReadId( );
  vec<longlong> prim2;
  vec<KmerPath> ext2;
  AllExtensions( jpaths[id2], scloud, local_ids, prim2, ext2 );

  // Orientations.
  bool rc1 = false;
  bool rc2 = false;
  if ( ins_db.IsValidType( insert_id ) ) rc2 = true;
  else {
    if ( ins_db.ReadLoc( insert_id ).Rc( ) ) rc1 = true;
    else rc2 = true;
  }
    
  // Try to walk all possible combinations
  for (int ii=0; ii<prim1.isize( ); ii++) {
    for (int jj=0; jj<prim2.isize( ); jj++) {
      int u1 = cloud_pathsdb[prim1[ii]].ReadId( );
      int u2 = cloud_pathsdb[prim2[jj]].ReadId( );
      if ( rc1 ) u1 = to_rc[u1];
      if ( rc2 ) u2 = to_rc[u2];
      
      vecKmerPath loc_c;
      vec<closure_t> loc_ct;
      AssembleCloud( max_klen, max_depth, verbose_log, u1, u2,
		     scloud, local_ids, ins_db, loc_c, loc_ct, log );
      
      closures.Append( loc_c );
      copy( loc_ct.begin( ), loc_ct.end( ), back_inserter( closure_types ) );
    }
  }
  
}

/**
 * TrimSeedingContigs
 */
void TrimSeedingContigs( const int u1,
			 const int u2,
			 const vecKmerPath &cloud_paths,
			 vecKmerPath &closures,
			 vec<closure_t> &closure_types )
{
  vecKmerPath new_closures;
  vec<closure_t> new_closure_types;
  new_closures.reserve( closures.size( ) );
  new_closure_types.reserve( closure_types.size( ) );
  
  // Nothing to do.
  if ( closures.size( ) < 1 ) return;

  // Loop over all closures.
  for (size_t ii=0; ii<closures.size( ); ii++) {
    const closure_t ctype = closure_types[ii];
    const KmerPath &ckmers = closures[ii];
    const KmerPath &u1kmers = cloud_paths[u1];
    const KmerPath &u2kmers = cloud_paths[u2];
    
    bool ltrim = ( ctype == FULL_CLOSURE || ctype == LEFT_CLOSURE );
    bool rtrim = ( ctype == FULL_CLOSURE || ctype == RIGHT_CLOSURE );

    int first_seg = 0;
    int last_seg = ckmers.NSegments( ) - 1;
    if ( ltrim ) first_seg += u1kmers.NSegments( ) - 1;
    if ( rtrim ) last_seg += - u2kmers.NSegments( ) + 1;

    // Check perfect matches.
    if ( ltrim ) {
      for (int jj=0; jj<first_seg; jj++)
	ForceAssert( ckmers.Segment( jj ) == u1kmers.Segment( jj ) );
      const KmerPathInterval &cseg = ckmers.Segment( first_seg );
      const KmerPathInterval &u1seg = u1kmers.Segment( first_seg );
      ForceAssert( cseg.Start( ) <= u1seg.Start( ) );
      ForceAssert( cseg.Stop( ) >= u1seg.Stop( ) );
    }

    if ( rtrim ) {
      const KmerPathInterval &cseg = ckmers.Segment( last_seg );
      const KmerPathInterval &u2seg = u2kmers.Segment( 0 );
      ForceAssert( cseg.Start( ) <= u2seg.Start( ) );
      ForceAssert( cseg.Stop( ) >= u2seg.Stop( ) );
      for (int jj=1; jj<u2kmers.NSegments( ); jj++)
	ForceAssert( ckmers.Segment( last_seg + jj ) == u2kmers.Segment( jj ) );
    }

    // Build trimmed KmerPath.
    KmerPath trimmed;
    for (int jj=first_seg; jj<=last_seg; jj++) {
      const KmerPathInterval &cseg = ckmers.Segment( jj );
      longlong start = cseg.Start( );
      longlong stop = cseg.Stop( );

      if ( ltrim && jj == first_seg ) {
	if ( stop == u1kmers.LastSegment( ).Stop( ) ) continue;
	start = u1kmers.Segment( jj ).Stop( ) + 1;
      }
      
      if ( rtrim && jj == last_seg ) {
	if ( start == u2kmers.Segment( 0 ).Start( ) ) continue;
	stop = u2kmers.Segment( 0 ).Start( ) - 1;
      }
      
      trimmed.AddSegment( start, stop );
    }

    // Not a proper closure.
    if ( trimmed.NSegments( ) < 1 ) continue;

    // Add closure to list.
    new_closures.push_back( trimmed );
    new_closure_types.push_back( ctype );
  }
  
  // Done.
  swap( closures, new_closures );
  swap( closure_types, new_closure_types );
  
}

/**
 * BvecPerfectOverlap
 *
 * NB: SmithWatBandedA does not find a long perfect align, even if it
 * exists, if band is bigger than the length of ( bases1 - offset )
 * (sante, 2013.04.08).
 */
int BvecPerfectOverlap( const int small_k,
			const int big_k,
			const bvec &bases1,
			const bvec &bases2 )
{
  align al;
  int err = 0;
  int offset = Max( 0, (int)bases1.size( ) - big_k );
  int band = Min( big_k, (int)bases1.size( ) - offset - small_k + 1 );
  SmithWatBandedA( bases1, bases2, offset, band, al, err );
  
  if ( al.Pos1( ) != (int)bases1.size( ) ) return 0;
  if ( ! al.Perfect( bases1, bases2 ) ) return 0;
  if ( al.Pos1( ) - al.pos1( ) < small_k ) return 0;

  return al.Pos1( ) - al.pos1( );
}

/**
 * MergePartialClosures
 */
void MergePartialClosures( const int small_k,
			   const int big_k,
			   CClosures &closures,
			   ostream *log )
{
  const vecbvec &bases_in = closures.Bases( );
  const vec<closure_t> &ctypes_in = closures.Ctypes( );

  CClosures closures_out;

  // Keep track of unused closures.
  vec<bool> used( bases_in.size( ), false );

  // First loop, which may leave out some partial closures.
  closures_out.Reserve( closures.Size( ) );
  for (int ii=0; ii<(int)bases_in.size( ); ii++) {
    if ( used[ii] ) continue;
    
    // An overlap or already full closure, add it to output.
    if ( ctypes_in[ii] == OVERLAP ) {
      closures_out.AddOverlap( closures.Overlap( ii ) );
      used[ii] = true;
      continue;
    }
    if ( ctypes_in[ii] == FULL_CLOSURE ) {
      closures_out.AddTrueClosure( closures.Bases( ii ), FULL_CLOSURE );
      used[ii] = true;
      continue;
    }

    // Skip this case for now.
    if ( ctypes_in[ii] == RIGHT_CLOSURE )
      continue;

    // A partial left closure: test against all partial right closures.
    ForceAssert( ctypes_in[ii] == LEFT_CLOSURE );
    for (int jj=0; jj<(int)bases_in.size( ); jj++) {
      if ( used[ii] ) break;
      if ( used[jj] ) continue;
      if ( ctypes_in[jj] != RIGHT_CLOSURE ) continue;
      
      if ( log )
	*log << "Testing partial left closure " << ii
	     << " (of size" << bases_in[ii].size( ) << ")\n";
      
      const bvec &bases1 = bases_in[ii];
      const bvec &bases2 = bases_in[jj];
      int overlap = BvecPerfectOverlap( small_k, big_k, bases1, bases2 );
      if ( overlap < small_k ) continue;
      
      // Merge partial closures.
      int al_pos1 = (int)bases_in[ii].size( ) - overlap;
      bvec frag( bases_in[ii], 0, al_pos1 );
      bvec merged = Cat( frag, bases_in[jj] );
      closures_out.AddTrueClosure( merged, FULL_CLOSURE );
      used[ii] = true;
      used[jj] = true;
      
      if ( log )
	*log << "Merged partial closures of sizes " << bases_in[ii].size( )
	     << " and " << bases_in[jj].size( )
	     << " into full closure of size " << merged.size( )
	     << "\n";
    }
    
  }

  // Take care of the still unused partial closures.
  for (int ii=0; ii<(int)bases_in.size( ); ii++)
    if ( ! used[ii] )
      closures_out.AddClosure( closures, ii );

  // Swap and done.
  swap( closures, closures_out );

}

/**
 * FillGap
 */
bool FillGap( const int small_k,
	      const int big_k,
	      const int max_klen,
	      const int max_depth,
	      const int max_dist,
	      const bool verbose_log,
	      const int gap_id,
	      const int u1,
	      const int u2,
	      const SCloud &scloud,
	      const CInsertsDB &ins_db,
	      const vecKmerPath &jpaths,
	      const KmerBaseBroker &kbb,
	      CClosures &closures,
	      ostream *log )
{
  closures.Clear( );
  
  const vecKmerPath &cloud_paths = scloud.paths_;
  const vec<tagged_rpint> &cloud_pathsdb = scloud.pathsdb_;
  const bvec u1_bases = kbb.Seq( cloud_paths[u1] );
  const bvec u2_bases = kbb.Seq( cloud_paths[u2] );

  // Log info.
  if ( log )
    *log << "FillGap( ) between " << u1 << " (" << u1_bases.size( ) << " bases)"
	 << " and " << u2 << " (" << u2_bases.size( ) << " bases)\n";
  
  // Test for direct overlap.  
  int overlap = BvecPerfectOverlap( small_k, big_k, u1_bases, u2_bases );
  if ( overlap > 0 ) {
    closures.AddOverlap( overlap );
    if ( log ) *log << "found perfect overlap of " << overlap << " bases\n";
  }

  // Localize global reads.
  vec<int> loc_ids;
  vec< pair<int,int> > tmp;
  LocalizedGlobalCloud( gap_id, u1, u2, max_dist, scloud, ins_db, jpaths, tmp );
  loc_ids.reserve( loc_ids.size( ) + tmp.size( ) );
  for (size_t jj=0; jj<tmp.size( ); jj++)
    loc_ids.push_back( tmp[jj].first );
  
  if ( log ) *log << "global cloud size: " << loc_ids.size( ) << "\n";
  
  // Fill gap using localized global reads.
  vecKmerPath kpaths;
  vec<closure_t> ktypes;
  AssembleCloud( max_klen, max_depth, verbose_log,
		 u1, u2, scloud, loc_ids, ins_db, kpaths, ktypes, log );
  
  // Remove seeding contigs from closures.
  TrimSeedingContigs( u1, u2, cloud_paths, kpaths, ktypes );
  
  // Try to merge partial closures.
  closures.Reserve( closures.Size( ) + kpaths.size( ) );
  for (int ii=0; ii<(int)kpaths.size( ); ii++)
    closures.AddTrueClosure( kbb.Seq( kpaths[ii] ), ktypes[ii] );
  
  ostream *mlog = verbose_log ? log : 0;
  MergePartialClosures( small_k, big_k, closures, mlog );
  
  // Done
  if ( log ) closures.PrintReport( *log );
  return ( closures.Size( ) > 0 );

}

/**
 * SelectClosures
 *
 * K: needed to define min gap size
 * gap_size: esimated gap size
 * gap_multiplier: needed to define max gap size
 * all_full: if true, select all full closures (within min-max range)
 */
void SelectClosures( const int K,
		     const int gap_size,
		     const double gap_multiplier,
		     const bool all_full,
		     CClosures &closures )
{
  // Heuristics.
  int min_gap = - K + 1;
  int max_gap = (int)( Max( 100.0, double( gap_size ) ) * gap_multiplier );

  // Nothing to select.
  if ( closures.Size( ) < 1 ) return;

  // The output (only full closures or overlaps, for now).
  CClosures out_closures;
  out_closures.Reserve( closures.Size( ) );
  for (size_t ii=0; ii<closures.Size( ); ii++) {

    // Not a full closure or overlap.
    if ( closures.IsLeft( ii ) || closures.IsRight( ii ) )
      continue;

    // Overlap.
    if ( closures.IsOverlap( ii ) ) {
      int over = closures.Overlap( ii );
      if ( -over >= min_gap ) 
	out_closures.AddOverlap( closures.Overlap( ii ) );
      continue;
    }

    // Full true closure.
    if ( closures.IsFull( ii ) ) {
      int closure_len = closures.Len( ii );
      if ( min_gap <= closure_len && closure_len <= max_gap )
	out_closures.AddClosure( closures, ii );
      continue;
    }

    // Should not get here.
    ForceAssert( 1 == 0 );
  }

  // If not all_full, just return the first closure.
  if ( ( ! all_full ) && out_closures.Size( ) > 1 )
    out_closures.Resize( 1 );

  // Swap and done.
  swap( closures, out_closures );
  
}

/**
 * SaveFilledAssembly
 */
void SaveFilledAssembly( const String &closures_head,
			 const String &assembly_in_head,
			 const String &assembly_out_head )
{
  const String closures_fastb_file = closures_head + ".fastb";
  const String closures_info_file = closures_head + ".info";
  const String in_contigs_fastb_file = assembly_in_head + ".contigs.fastb";
  const String in_superb_file = assembly_in_head + ".superb";
  
  // Load closures.
  vecbvec cl_bases( closures_fastb_file );
  int n_closures = cl_bases.size( );
  
  vec< quad<int,int,int,char> > cl_info;
  cl_info.reserve( n_closures );
  int sid = -1;
  int gid = -1;
  int over = 0;
  char cinfo = 'x';
  ifstream cl_in( closures_info_file.c_str( ) );
  while ( cl_in ) {
    cl_in >> sid >> gid >> over >> cinfo;
    if ( ! cl_in ) break;
    cl_info.push_back( quad<int,int,int,char>( sid, gid, over, cinfo ) );
  }
  cl_in.close( );
  sort( cl_info.begin( ), cl_info.end( ) );
  
  // Load assembly.
  vecbvec in_cbases( in_contigs_fastb_file );
  vec<superb> in_supers;
  ReadSuperbs( in_superb_file, in_supers );
  int n_supers = in_supers.size( );
  
  // Output contigs and supers.
  vecbvec out_cbases;
  vec<superb> out_supers;
  out_cbases.reserve( in_cbases.size( ) );
  out_supers.reserve( n_supers );
  
  // Loop over all supers.
  vec< quad<int,int,int,char> >::iterator qit;
  vec< quad<int,int,int,char> >::iterator qit2;
  for (int super_id=0; super_id<n_supers; super_id++) {
    const superb &insup = in_supers[super_id];
    const int n_intigs = insup.Ntigs( );

    // First contig in super.
    const bvec &first_tig = in_cbases[ insup.Tig( 0 ) ];
    int first_len = first_tig.size( );
    int first_id = out_cbases.size( );
    
    superb new_super;
    out_cbases.push_back( first_tig );
    new_super.PlaceFirstTig( first_id, first_len );
    
    // Loop over all contigs in insup (each iter: add prev gap, then contig).
    for (int tpos=1; tpos<n_intigs; tpos++) {
      bvec tig_bases = in_cbases[insup.Tig( tpos )];   // local copy
      int tig_len = tig_bases.size( );
      int tig_id = out_cbases.size( );
      int gap_size = insup.Gap( tpos - 1 );
      int gap_dev = insup.Dev( tpos - 1 );
      
      // Find closures for gap before tpos.
      quad<int,int,int,char> qgap( super_id, tpos-1, 0, '0' );
      qit = lower_bound( cl_info.begin( ), cl_info.end( ), qgap );
      vec< quad<int,int,int,char> > closures;
      vec<int> closure_ids;
      for ( qit2=qit; qit2!=cl_info.end( ); qit2++) {
	if ( qit2->first != super_id || qit2->second != tpos-1 ) break;
	closures.push_back( *qit2 );
	closure_ids.push_back( distance( cl_info.begin( ), qit2 ) );
      }
      
      ForceAssert( closures.size( ) < 3 );
      if ( closures.size( ) == 2 ) {
	char cl_type0 = closures[0].fourth;
	char cl_type1 = closures[1].fourth;
	ForceAssert( cl_type0 == 'L' && cl_type1 == 'R' );
      }
      
      // No closures found.
      if ( closures.size( ) == 0 ) {
	out_cbases.push_back( tig_bases );
	new_super.AppendTig( tig_id, tig_len, gap_size, gap_dev );
	continue;
      }

      // Full closure or overlap found.
      if ( closures.size( ) == 1 ) {
	int cl_id = closure_ids[0];
	int cl_size = cl_bases[cl_id].size( );
	
	bvec addon;
	if ( closures[0].fourth == 'O' ) {
	  ForceAssertEq( cl_size, 0 );
	  int overlap = closures[0].third;
	  addon = bvec( tig_bases, overlap, tig_len - overlap );
	}
	else if ( closures[0].fourth == 'F' ) {
	  ForceAssertGt( cl_size, 0 );
	  addon = Cat( cl_bases[cl_id], tig_bases );
	}
	else {
	  // Can be L or R.
	}

	if ( addon.size( ) > 0 ) {
	  bvec &out_last = out_cbases[out_cbases.size( )-1];
	  out_last = Cat( out_last, addon );
	  continue;
	}
      }

      // Partial closures found.
      for (int tt=0; tt<closures.isize( ); tt++) {
	int cl_id = closure_ids[tt];
	int cl_size = cl_bases[cl_id].size( );
	ForceAssert( closures[tt].fourth == 'L' || closures[tt].fourth == 'R' );
	ForceAssertGt( cl_size, 0 );
	
	if ( closures[tt].fourth == 'L' ) {
	  bvec &out_last = out_cbases[out_cbases.size( )-1];
	  out_last = Cat( out_last, cl_bases[cl_id] );
	}
	else {
	  tig_bases = Cat( cl_bases[cl_id], tig_bases );
	  tig_len = tig_bases.size( );
	}

	gap_size += -cl_size;
      }

      out_cbases.push_back( tig_bases );
      new_super.AppendTig( tig_id, tig_len, gap_size, gap_dev );
      
    } // loop over all contigs in insup

    // Adjust lengths of contigs in new_super, and append it to out_supers.
    for (int ii=0; ii<new_super.Ntigs( ); ii++)
      new_super.SetLen( ii, (int)out_cbases[new_super.Tig( ii )].size( ) );

    out_supers.push_back( new_super );

  } // loop over all supers

  // Convert to fastavector.
  vec<fastavector> out_fastavec;
  out_fastavec.reserve( out_cbases.size( ) );
  for (size_t ii=0; ii<out_cbases.size( ); ii++)
    out_fastavec.push_back( fastavector( out_cbases[ii] ) );

  // And save.
  SaveScaffoldAssembly( assembly_out_head, out_supers, out_fastavec, 0, 1 );
  
}

