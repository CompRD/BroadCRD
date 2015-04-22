///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__ASSISTED__CLOUD_WALKER__H
#define PATHS__ASSISTED__CLOUD_WALKER__H

#include "Basevector.h"
#include "SeqInterval.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"
#include "paths/KmerPathInterval.h"
#include "paths/assisted/CClosures.h"
#include "paths/assisted/CWalk.h"
#include "paths/assisted/CInsertsDB.h"
#include "paths/assisted/KmerPathPerfectOverlaps.h"

/**
 * struct SCloud
 *
 * Container for the data structure needed to represent a global cloud.
 */
struct SCloud {

  SCloud( const vec<int> &to_rc, 
	  const vecKmerPath &paths,
	  const vec<tagged_rpint> &pathsdb ) :
    to_rc_ ( to_rc ), paths_ ( paths ), pathsdb_( pathsdb ) { }

  const vec<int> &to_rc_;
  const vecKmerPath &paths_;
  const vec<tagged_rpint> &pathsdb_;
  
};

/**
 * OverlappingBuddies
 *
 * Find all sequences in the glocbal cloud that overlap the given kmer
 * path.
 */
void OverlappingBuddies( const SCloud &scloud,
			 const KmerPath &kpath,
			 vec<int> &buddies );

/**
 * SecondaryCloud
 *
 * Append secondary cloud of sequences aligning existing cloud. It
 * returns the number of sequences added to ids2dist.
 */
int SecondaryCloud( const int level,
		    const SCloud &scloud,
		    vec< pair<int,int> > &ids2dist );

/**
 * LocalizedGlobalCloud
 *
 * Build the list of ids of sequences from the global cloud local to
 * the selected gap. It first find the cloud of jump mates local to
 * gap_id, and then the ids of all the sequences in the global cloud
 * that overlap the localized jump mates.
 *
 * Output consists of pairs ids2dist, where first contains the id of
 * the sequence in scloud.paths_, and dist the distance from a seeding
 * jump. 0 means the sequence aligns directly to a jump read; i that
 * the sequence aligns directly to a sequence with dist (i-1), and to
 * no sequence with dist < (i-1).
 */
void LocalizedGlobalCloud( const int gap_id,
			   const int gcid1,
			   const int gcid2,
			   const int MAX_DIST,
			   const SCloud &scloud,
			   const CInsertsDB &ins_db,
			   const vecKmerPath &jpaths,
			   vec< pair<int,int> > &ids2dist );

/**
 * AllExtensions
 *
 * Find the locs of all the cloud sequences extending maximally the
 * given read (with perfect overlaps with the read), from a localized
 * set of sequences from a global cloud:
 *
 *          the read
 *       --->
 *    u1  ------------x-
 *    u2   -----------x----
 *    u3   -------*----------
 *
 * In the example above, the code will return u2 and u3, but not u1,
 * since u2 is a longer extension than u1 (and u1 matches u2). Note
 * that u1, u2, and u3 are required to align perfectly the read (and
 * all of them belong to a local cloud, as specified by select).
 *
 * NB: the code returns both the locs of the longest extensions (as
 *     indexes in the cloud_pathsdb), and the actual extensions (as
 *     KmerPaths).
 *
 * NB: the code relies on the assumption that sequences in the global
 *     cloud appear both as fw and rc.
 */
void AllExtensions( const KmerPath &read,
		    const SCloud &scloud,
		    const vec<int> &select,
		    vec<longlong> &extensions,
		    vec<KmerPath> &extended_paths );

/**
 * CloudWalk
 *
 * Iteratively walk from one sequence to the other. This is the
 * recursive core of the cloud walking algorithm, and it is called by
 * CloudWalkWrapper.
 */
void CloudWalk( const int max_klen,
		const int max_depth,
		const int target_id,
		const SCloud &scloud,
		const vec<int> &select,		
		int depth,
		vec<CWalk> &walks,
		ostream *log = 0 );

/**
 * CloudWalkWrapper
 *
 * Entry point for the core recursive function (CloudWalk). This is an
 * utility function called by AssembledCloud.
 */
void CloudWalkWrapper( const int u1,
		       const int u2,
		       const int max_klen,
		       const int max_depth,
		       const bool verbose_log,
		       const SCloud &scloud,
		       const vec<int> &local_ids,
		       vec<CWalk> &walks,
		       ostream *log = 0 );

/**
 * AssembleCloud
 *
 * Assemble local cloud. Output is stored in closures (partial
 * closures allowed).
 *
 * There are two versions of AssembleCloud, the difference being that
 * in the second one we close an insert, in the first we close the gap
 * between two contigs (not necessarily adjacent) in the same super
 * (note that the second version calls the first version multiple
 * times. I should probably get rid of the second version altogether).
 *
 * NB: closures will always contain the kmers of the seeding contigs,
 *    or of the contigs containing the seeding insert ends. If only
 *    partial closures are found, then each partial closure will
 *    contain the kmers of one of the contigs.
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
		    ostream *log = 0 );

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
		    ostream *log = 0 );

/**
 * TrimSeedingContigs
 *
 * Remove the seeding contigs from the given closures. Note that this
 * may reduce the size of closures (and closure_types), possibly to
 * zero.
 */
void TrimSeedingContigs( const int u1,
			 const int u2,
			 const vecKmerPath &cloud_paths,
			 vecKmerPath &closures,
			 vec<closure_t> &closure_types );

/**
 * BvecPerfectOverlap
 *
 * Decide if the two bvecs overlap perfectly by at least small_k
 * bases. It uses big_k to define the offset, and it returns the
 * amount of the overlap (or 0 if not found).
 *
 * NB: all aligns in which bases2 does not extend bases1 are
 * discarded!
 */
int BvecPerfectOverlap( const int small_k,
			const int big_k,
			const bvec &bases1,
			const bvec &bases2 );
  
/**
 * MergePartialClosures
 *
 * Try to merge partial closures, by allowing perfect matches >=
 * small_k.
 */
void MergePartialClosures( const int small_k,
			   const int big_k,
			   CClosures &closures,
			   ostream *log = 0 );

/**
 * FillGap
 *
 * Fill gap between contigs u1 and u2 (where win describes the gap
 * interval in the super it belongs to).  Output consists of closures
 * (partial or full), and/or of an "overlap" amount (ie, amount of
 * overlap between u1 and u2).  It returns false if no closure is
 * found.
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
	      ostream *log = 0 );

/**
 * SelectClosures
 *
 * Select winner closure(s) from the output of FillGap.
 */
void SelectClosures( const int K,
		     const int gap_size,
		     const double gap_multiplier,
		     const bool all_full,
		     CClosures &closures );

/**
 * SaveFilledAssembly
 *
 * Load closures, apply closures and save assembly.
 */
void SaveFilledAssembly( const String &closures_head,
			 const String &assembly_in_head,
			 const String &assembly_out_head );

#endif
