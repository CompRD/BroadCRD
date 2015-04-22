/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * JumpAligner: this is a highly modified adaptation of
 * FirstLookupFinder which is specific to finding the best alignments
 * for jumping reads, which may include a junction crossing (crossing
 * over to a distance sequence).  We want to find the most plausible
 * home(s).  We pass in reads, quals, unipaths, adjacency graph, and
 * we return the best matching sequence(s) for the read, even those
 * which cross unipath boundaries through the graph.  --bruce 4 Dec 09
 *
 */
#include "VecUtilities.h"
#include "lookup/JumpAligner.h"
#include "system/WorklistN.h"

#define IMITATE_BUG 1
#define DUMP_BASES 0
#if DUMP_BASES
#include <iostream>
template <class Itr>
void dumpBases( Itr begin, Itr const& end )
{ while ( begin != end )
  { std::cout << Base::val2Char(*begin); ++begin; }
  std::cout << std::endl; }
#endif

#define is_in(x,v) (find(v.begin(), v.end(), x) != v.end())

namespace
{

class Processor
{
public:
    Processor( JumpAligner const& flf, vecbvec const& queries, vecqvec const& quals,
               vec<jump_align>& alignments, LockedData& lock,
               u_int chunkSize )
      : mFLF(flf), mQueries(queries), mQuals(quals), mAlignments(alignments), mLock(lock),
	mChunkSize(chunkSize)
    {}

    Processor( Processor const& that )
      : mFLF(that.mFLF), mQueries(that.mQueries), mQuals(that.mQuals), mAlignments(that.mAlignments),
	mLock(that.mLock), mChunkSize(that.mChunkSize)
    {}

    // copy-assignment prohibited by reference members, compiler-supplied destructor is OK

    void operator()( u_int batchNo ) const
    {
        u_int queryID = batchNo*mChunkSize;
        u_int end = std::min(queryID+mChunkSize,
                                    static_cast<u_int>(mQueries.size()));
	std::list<jump_align> results;
        for ( ; queryID < end; ++queryID )
        {
	  //vec<jump_align> aligns;
	  mFLF.getAlignments(mQueries[queryID],mQuals[queryID], queryID, results );
	  //results.insert(results.end(),aligns.begin(),aligns.end());
        }
        if ( results.size() )
        {
            Locker lock(mLock);
            mAlignments.insert(mAlignments.end(),results.begin(),results.end());
        }
    }

private:
    JumpAligner const& mFLF;
    vecbvec const& mQueries;
    vecqvec const& mQuals;
    vec<jump_align>& mAlignments;
    LockedData& mLock;
    u_int mChunkSize;
};

}

// trim path to n bases
void trim_path(jump_align_path & jap, u_int new_length) {
  for (u_int s = 0, n = new_length; s < jap.size(); ++s) {
    if (n <= jap[s].length) {
      jap[s].length = n;
      jap.resize(s+1);
      break;
    } else n -= jap[s].length;
  }
}


// see if these are really the same path
bool same_path(jump_align_path const& p1, jump_align_path const& p2, u_int mtl) {
  u_int psize = p1.size();
  if (psize != p2.size()) return false;
  for (u_int s = 0; s < psize; ++s) {
    if (p1[s].contig != p2[s].contig
	|| p1[s].offset != p2[s].offset
	|| p1[s].length != p2[s].length) {
      return false;
    }
  }
  return true;
}


void JumpAligner::getAlignments( bvec const& query,
				 qvec const& qual,
				 u_int queryID,
				 std::list<jump_align> & result ) const
{
#if DUMP_BASES
std::cout << "Query: " << queryID << std::endl;
#endif
    u_int qsize = query.size();
    u_int lookupK = mLookupTab.getK();
    u_int targetLength = max(qsize, mFilter.target_length);
    bool debug = is_in(queryID, mFilter.debug_reads);

    if (query.size() >= lookupK)
    {
      vec<jump_align_path> targets;

      if ( mFilter.orientation != JumpAlignFilter::FW_ONLY ) {
	// Get RC candidates
	u_int kStart = qsize - lookupK;
	u_int kmer = mLookupTab.getKmer(query.RCBegin(kStart));
	LocationVec const& locations = mLookupTab.getLocations(kmer);
	extendLocations(locations, targetLength, True, targets, debug);
      }
      if ( mFilter.orientation != JumpAlignFilter::RC_ONLY ) {
	// Get FW candidates
        u_int kmer = mLookupTab.getKmer(query.Begin());
        LocationVec const& locations = mLookupTab.getLocations(kmer);
	extendLocations(locations, targetLength, False, targets, debug);
      }
      u_int ntargets = targets.size();
      if (ntargets) {
	u_int mtl = maxTrustedLength(query, targets, debug);
	if (debug) PRINT3(queryID, mtl, ntargets);
	if (mtl >= mFilter.min_match) {
	  // score and sort targets best to worst
	  vec<u_int> scores(ntargets, (u_int) 1000000);
	  scoreAlignments(query, qual, mtl, targets, scores, debug);
	  if (ntargets > 1) SortSync(scores, targets);
	  int bestscore = (int) scores[0];
	  if (mFilter.score_max < 0 || bestscore <= mFilter.score_max) {
	    u_int good = 0;
	    u_int max_to_return = ntargets;
	    if (mFilter.max_placements > 0 && mFilter.max_placements < max_to_return)
	      max_to_return = mFilter.max_placements;
	    //u_int score_limit = bestscore + mFilter.score_delta;
	    // if there is a perfect score, don't return imperfect alternatives
	    u_int score_limit = (bestscore==0) ? 0 : (bestscore + mFilter.score_delta);
	    for (u_int i = 0; i < ntargets && good < max_to_return && scores[i] <= score_limit; ++i) {
	      jump_align_path &t = targets[i];
	      jump_align_path &last_path = i ? targets[i-1] : targets[i];
	      jump_align ja;
	      ja.query_id = queryID;
	      ja.query_length = qsize;
	      trim_path(t, mtl);
	      if (i > 0 && same_path(t, last_path, mtl)) {
		if (debug) {
		  PRINT4(scores[i], t[0].contig, last_path[0].contig, "dup");
		}
		continue;
	      } else {
		last_path = t;
		ja.target_path = t;
		ja.rc = (t[0].offset < 0);
		ja.trusted_length = mtl;
		if (ja.rc) ja.unReverse(mContigs);
		ja.extended = (ja.target_path.size() > 1);
		result.push_back(ja);
		++good;
		if (debug) {
		  PRINT3(scores[i], t[0].contig, t[0].offset);
		}
	      }
	    }
	  }
	}
      }
    }
}

void JumpAligner::getAllAlignments( vecbvec const& queries,
					  vecqvec const& quals,
                                          vec<jump_align>& alignments,
                                          u_int nThreads ) const
{
    // Reserve memory for alignments.  We make a heuristic estimate for the
    // upper bound of the alignment frequency per read.
    static const double ALIGNS_PER_READ = 1.25;
    alignments.reserve(queries.size() * ALIGNS_PER_READ);

    LockedData lock;
    Processor processor(*this, queries, quals, alignments, lock, CHUNK_SIZE);
    u_int nnn = queries.size();
    parallelFor(0u,(nnn+CHUNK_SIZE-1)/CHUNK_SIZE,processor,nThreads);
}




void
JumpAligner::extendThruGraph(u_int c,
			     u_int offset,
			     u_int length,
			     const bool rc,
			     vec<jump_align_path>& paths
			     ) const
{
  vec<int> nodes;
  u_int nbases;
  const bvec& contig = mContigs[c];
  u_int csize = contig.size();
  jump_align_path path;
  jump_align_segment segment;

  ForceAssertLt(offset, csize);

  nbases = min(length, csize-offset);
  segment.contig = c;
  segment.offset = rc ? -(offset+1) : offset;
  segment.length = nbases;
  if (rc) {
    nodes = mGraph.To(c);
  } else {
    nodes = mGraph.From(c);
  }
  if (nbases == length) {
    path.push_back(segment);
    paths.push_back(path);
  } else if (mFilter.cross_unipath_boundaries) {
    //PRINT4(csize, offset, nbases, length);
    // recurse thru the graph, skipping over the first K-1 bases of the next link
    for (u_int n = 0; n < nodes.size(); ++n) {
      vec<jump_align_path> tpaths;
      extendThruGraph(nodes[n], mK-1, length-nbases, rc, tpaths);
      for (u_int t = 0; t < tpaths.size(); ++t) {
	tpaths[t].push_front(segment);
	paths.push_back(tpaths[t]);
	if (paths.size() >= mFilter.max_kmer_freq) break;
      }
    }
  }
  //PRINT(results.size());
}



void 
JumpAligner::extendLocations(LocationVec const& locations,
			     u_int length,
			     const bool rc,
			     vec<jump_align_path>& paths,
			     const bool debug) const
{
  u_int nlocs = locations.size();
  if (nlocs > mFilter.max_kmer_freq) return;
  paths.reserve(floor(nlocs*1.2)); // make room for 
  for (u_int i = 0; i < nlocs; ++i) {
    Location loc = locations[i];
    u_int c = loc.getContig();
    bvec const& contig = mContigs[c];
    // offset into contig, from direction we are moving, i.e.,
    // how many bases to skip over (rc is right to left)
    u_int offset = rc ? (contig.size() - (loc.getOffset() + mLookupTab.getK())) : loc.getOffset();
    if (debug) PRINT4(i, c, contig.size(), offset)
    if (contig.size() - offset < mK) continue;
    extendThruGraph(c, offset, length, rc, paths);
    if (debug) PRINT(paths.size());
    if (paths.size() > mFilter.max_kmer_freq) {
      paths.clear();
      return;
    }
  }
}


u_int 
JumpAligner::maxTrustedLength(bvec const& query,
			      vec<jump_align_path> const& targets,
			      const bool debug) const
{
  u_int mtl = 0;
  u_int qsize = query.size();
  u_int nhood_size = mFilter.mismatch_neighborhood;
  //query.Print(cout);
  for (u_int t = 0; t < targets.size(); ++t) {
    jump_align_path_iterator japi(targets[t], mContigs, 0);
    u_int trustedLength = qsize; // until proven otherwise
#define NHOOD_BITS (8*sizeof(int)) // assuming 8-bit bytes
    ForceAssertLe(nhood_size, NHOOD_BITS);
    bitset<NHOOD_BITS> nhood;
    u_int mismatches = 0;
    for (u_int i = 0, n = 0; i < qsize; ++i, ++japi) {
      //u_int qi = query[i], ji = *japi;
      //PRINT3(i, qi, ji);
      int mismatch = (query[i] == *japi) ? 0 : 1;
      //if (mismatch) cout << "*";
      //else cout << ".";
      // remember old value in this slot before replacing it
      int was = nhood[n];
      nhood.set(n, mismatch);
      // update running total of mismatches in nhood
      mismatches += mismatch - was;
      // advance to next slot
      if (++n >= nhood_size) n = 0;
      // see if we exceed our limit in nhood
      if (mismatches >= mFilter.mismatch_threshhold) {
	trustedLength = i + 1;
	// go back before first mismatch within neighborhood
	// n is pointing at least recent slot
	for (u_int m = 0; m < nhood_size; ++m) {
	  u_int ni = n + m;
	  if (ni >= nhood_size) ni -= nhood_size;
	  if (nhood[ni]) {
	    trustedLength -= nhood_size - m;
	    break;
	  }
	}
	trustedLength -= min(trustedLength, mFilter.mismatch_backoff);
	break;
      }
    }
    if (debug) PRINT4(targets[t][0].contig, targets[t][0].offset, targets[t][0].length, trustedLength);
    //cout << endl;
    if (trustedLength > mtl)
      mtl = trustedLength;
    if (mtl == qsize) break;	// can't do better than this!
  }
  return mtl;
}


// 	scoreAlignments(query, qual, queryID, mtl, targets); Given a
// maximum trusted length from all candidate alignments, compute a
// score for each over that interval.  Score is just sum of quals at
// mismatch positions.  sort the targets best to worst.

void JumpAligner::scoreAlignments(bvec const& query,
				  qvec const& qual,
				  const u_int maxTrustedLength,
				  vec<jump_align_path>& targets,
				  vec<u_int>& scores,
				  const bool debug
			 ) const
{
  u_int qsize = query.size();
  u_int perfects = 0;
  ForceAssertLe(maxTrustedLength, qsize);
  for (u_int t = 0; t < targets.size(); ++t) {
    jump_align_path_iterator japi(targets[t], mContigs, 0);
    u_int score = 0;
    u_int mismatches = 0;
    u_int segn = 0, seglen = debug ? targets[t][segn].length : (u_int)-1;
    for (u_int i = 0; i < maxTrustedLength; ++i, ++japi) {
      if (query[i] != *japi) {
	score += qual[i];
	++mismatches;
	if (debug) cout << "*";
      } else if (debug) cout << ".";
      if (debug && seglen == i+1) {
	cout << '|';
	seglen += targets[t][++segn].length;
      }
    }
    scores[t] = score;
    if (debug) {
      cout << endl;
      PRINT3(score,mismatches,maxTrustedLength);
    }
    if (score == 0 && mFilter.max_placements && ++perfects >= mFilter.max_placements)
      break;
  }
}

