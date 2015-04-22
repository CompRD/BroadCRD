////////////////////////////////////////////////////////////////////////////
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
#ifndef LOOKUP_JUMPALIGNER_H_
#define LOOKUP_JUMPALIGNER_H_

#include "lookup/LookupTab.h"
#include "lookup/LookAlign.h"
#include "graph/Digraph.h"
#include "Basevector.h"
#include "Vec.h"
#include <list>
#include <algorithm>
#include <cstddef>

// jump_align is similar to other alignement structures, but it also
// includes the actual bvec of bases with the best match, since they
// may be from branching unipaths.  We also include varous stats from
// the alignment process.


struct jump_align_segment {
  u_int contig, length;
  int offset;			/* negative means rc (from other end) */
};

typedef vec<jump_align_segment> jump_align_path; // contig, offset, length

class jump_align {
public:
  unsigned int query_id, query_length, score, trusted_length;
  jump_align_path target_path;
  bool rc, extended;
  
  // reverse the path through the segments, updating offsets/lengths
  void unReverse(const vecbasevector & target) {
    ForceAssertEq(this->rc, True);
    unsigned int path_length = target_path.size();
    //bool debug = (path_length > 1);
    bool debug = false;
    if (debug) {
      PRINT(query_id);
      for (unsigned int i = 0; i < path_length; ++i) {
	jump_align_segment& s = target_path[i];
	unsigned int csize = target[s.contig].size();
	PRINT4(s.contig, s.offset, s.length, csize);
      }
      PRINT(path_length);
    }
    for (unsigned int i = 0, length = trusted_length, next_offset = 0; i < path_length; ++i) {
      unsigned int j = path_length - (i+1);
      if (i < j) swap(target_path[i], target_path[j]);
      jump_align_segment& s = target_path[i];
      unsigned int csize = target[s.contig].size();
      if (i == 0) {
	// figure out what the K-1 step is across unipaths
	// by looking at the old last seg offset (now first)
	next_offset = -s.offset - 1;
	// compute next offset from other end
	s.offset = csize - (s.length - s.offset) + 1;
      } else {
	s.offset = next_offset;
      }
      u_int n = csize - s.offset;
      if (length > n) s.length = n;
      else s.length = length;
      length -= s.length;
      if (debug) PRINT4(s.contig, s.offset, s.length, csize);
      // make sure we use the whole path
      if (i == path_length-1) {
	ForceAssertEq(length, 0u);
      } else if (length == 0) ForceAssertEq(i, path_length-1);
    }
    this->rc = False;
  }
};


class jump_align_path_iterator
: public std::iterator<std::input_iterator_tag,unsigned char,
  std::ptrdiff_t, void, unsigned char>
{
 public:
 jump_align_path_iterator( jump_align_path const& jap, vecbvec const& contigs, unsigned long pos = 0 )
   : mJAP(jap), mContigs(contigs), mPathLength(jap.size()), mRC(jap.size() > 0 && jap[0].offset < 0) {
    set_pos(pos);
  }

  // compiler-supplied copying and destructor are OK

  jump_align_path_iterator& operator++() {
    advance();
    return *this;
  }

  jump_align_path_iterator operator++( int ) {
    jump_align_path_iterator tmp(*this);
    advance();
    return tmp;
  }

  unsigned char operator*() const {
    u_int base;
    if (mRC) {
      base = BaseVec::Complement((*mSegment)[mSegmentBase - mBaseIndex]); /* right to left */
    } else {
      base = (*mSegment)[mSegmentBase + mBaseIndex]; /* left to right */
    }
    return base;
  }

  friend bool operator==( jump_align_path_iterator const& fi1, jump_align_path_iterator const& fi2 ) {
    return &fi1.mJAP == &fi2.mJAP &&
    fi1.mSegIndex == fi2.mSegIndex &&
    fi1.mBaseIndex == fi2.mBaseIndex;
  }


  friend bool operator!=( jump_align_path_iterator const& fi1, jump_align_path_iterator const& fi2 )
  {
    return !(fi1==fi2);
  }

  

 private:
  void adjust_segment() {
    const jump_align_segment &s = mJAP[mSegIndex];
    mSegment = &mContigs[s.contig];
    mSegmentBase = mRC ? ((*mSegment).size() + s.offset) : s.offset;
  }

  void advance(u_int n = 1) {
    mBaseIndex += n;
    const jump_align_segment &s = mJAP[mSegIndex];
    while (mBaseIndex >= s.length) {
      mBaseIndex -= s.length;
      ++mSegIndex;
      // if we've hit the end, make sure we don't go further
      if (mSegIndex == mPathLength) { ForceAssertEq(mBaseIndex, 0u); }
      else adjust_segment();
    }
  }

  void set_pos(size_t pos) {
    mSegIndex = 0;
    mBaseIndex = 0;
    adjust_segment();
    if (pos) advance(pos);
  }

  const jump_align_path& mJAP;
  const vecbvec& mContigs;
  u_int mPathLength;
  vecbvec::size_type mSegIndex; /* index into jump_align_path */
  bvec::size_type mBaseIndex;	 /* index int path segment */
  const bvec *mSegment;
  bvec::size_type mSegmentBase; /* start of segment within contig */
  bool mRC;
};


// JumpAlignFilter
// A struct to hold filtering options for a call to FirstLookup.  These
// filters are all applied early on in the FirstLookup algorithm and serve to
// reduce runtime.  The default behavior is to apply no filters.
struct JumpAlignFilter {

  // Constructor: Sets default values (i.e., no filtering).
  // [Actually, not true...the new scoring parameters do some filtering...
  // either need to clean up documentation or defaults.  --bruce]
  JumpAlignFilter( ) {
    orientation = ALL;
    max_kmer_freq = 0U;
    min_size = 0U;
    mismatch_threshhold = 3;
    mismatch_neighborhood = 8;
    mismatch_backoff = 3;
    max_error_rate = 1.f;
    min_match = 0;
    score_max = -1;
    score_delta = 0;
    target_length = 0;
    max_placements = 0;
    cross_unipath_boundaries = False;
    vec<int> debug_reads;
  }

  // Filtering options, listed in the order in which they are applied.

  // Only look for alignments where the orientation between query and target is as follows:
  enum { ALL, FW_ONLY, RC_ONLY } orientation;

  // Ignore reads whose leading kmers (using the lookup table's K) appear too
  // often in the lookup table.  This criterion is applied separately to the forward
  // and reverse directions before applying any other criteria.  A zero means that this
  // criterion is not applied.
  unsigned int max_kmer_freq;

  // Only extend from kmers where the potential alignment length is at least this great.
  // This is a criterion on how close to the end (beginning for reverse) of the
  // contig we landed, as well as the query length.
  unsigned int min_size;

  // These three parameters determine when to quit matching; if it sees mismatch_threshhold
  // incorrect bases in the most recent mismatch_neighborhood, it will then cut off the
  // match at "backoff" bases before the first bad base within the neighborhood.  This is deisgned
  // to identify jumping read junction crossings.
  unsigned int mismatch_threshhold;
  unsigned int mismatch_neighborhood;
  unsigned int mismatch_backoff;

  // If set, require this many bases to align (after cutoff and backoff described by the above
  // parameters) before accepting the alignment.
  unsigned int min_match;

  // Max error rate.
  float max_error_rate;

  // Sum of mismatched base qualities much not exceed this.  If negative, don't
  // check.
  int score_max;

  // Score delta: include all aligns with a score within this delta of
  // the best score (i.e., if this value is 0, only those alignment(s)
  // tied for the best score will be returned).  If positive, allows
  // for returning alignments with scores that close to the
  // best.
  unsigned int score_delta;

  // Extend results to this length (if greater than query length)
  unsigned int target_length;

  // If true, will look for read placements which cross unipath boundaries.
  Bool cross_unipath_boundaries;

  // max_placements: if greater than zero, don't bother returning more
  // than this number of placements even if they are close in score.
  unsigned int max_placements;

  // Not really a filtering option, but print debugging information about these reads, if set.
  vec<int> debug_reads;
};

class JumpAligner
{
public:
    /// constructor does not copy its args.  they must live as long as the finder.
 JumpAligner( JumpAlignFilter const& filter, LookupTab const& lookupTab, vecbvec const& contigs, digraph& graph, unsigned int K )
   : mFilter(filter), mLookupTab(lookupTab), mContigs(contigs), mGraph(graph), mK(K)
    {}

    // compiler-supplied copy construction and destructor are fine.
    // class is not assignable due to references.

    /// Returns a list of "first lookup" alignments for a specified query.
    /// These are the least-mismatched alignments extended from an initial perfectly matched kmer
    /// for the query or its reverse complement, and subject to the filtering criteria.
  void getAlignments( bvec const& query, qvec const& qual, unsigned int queryID, std::list<jump_align> & result ) const;

    /// The alignments vector is populated with all alignments for all the queries.
    /// The result will be identical to adding the results of calling getAlignments on each query in turn to
    /// the output vector, but the operation is actually done in parallel, so the results won't be in
    /// exactly the same order as if it had been done that way.
  void getAllAlignments( vecbvec const& queries, vecqvec const& quals, vec<jump_align>& alignments, unsigned int nThreads = 8 ) const;

private:
  void extendThruGraph(unsigned int node, unsigned int offset, unsigned int length, const bool rc, vec<jump_align_path>& paths
 ) const;
  void extendLocations(LocationVec const& locations, unsigned int length, const bool rc, vec<jump_align_path>& results, const bool debug) const;
  unsigned int maxTrustedLength(bvec const& query, vec<jump_align_path> const& targets, const bool debug) const;
  void scoreAlignments(bvec const& query, qvec const& qual, unsigned int matchLength, vec<jump_align_path>& alignments, vec<u_int>& scores, const bool debug) const;
  bool sortAlignments(const jump_align& a1, const jump_align& a2) const;

    JumpAlignFilter const& mFilter;
    LookupTab const& mLookupTab;
    vecbvec const& mContigs;
    digraph const& mGraph;
    // This is the passed in "big" K to determine unibase overlaps,
    // not the lookup table K
    unsigned int mK;
    static unsigned int const CHUNK_SIZE = 10000; // number of reads to process with one trip to the worklist
};

#endif /* LOOKUP_JUMPALIGNER_H_ */
