#include "pairwise_aligners/KmerAligner.h"

void FindPossibleOrientedAlignments(const int KMER_SIZE, 
                                    const basevector & target, 
                                    const vecbasevector & queries,
                                    vec<PerfectMatch> & aligns,
                                    const vec<int> & queryIds) {
  vecbasevector temp;
  temp.push_back(target);
  vec<int> targetIds(1,0);
  FindPossibleOrientedAlignments(KMER_SIZE, temp, queries, aligns,
                                 targetIds, queryIds);
}

void FindPossibleOrientedAlignments(const int KMER_SIZE, 
                                    const vecbasevector & target, 
                                    const vecbasevector & queries,
                                    vec<PerfectMatch> & aligns,
                                    const vec<int> & targetIds,
                                    const vec<int> & queryIds) { 
  switch (KMER_SIZE) {
    case 8: FindPossibleOrientedAlignments<8>
              (target,queries, aligns, targetIds, queryIds); break;
    case 12: FindPossibleOrientedAlignments<12>
               (target,queries, aligns, targetIds, queryIds); break;
    case 16: FindPossibleOrientedAlignments<16>
               (target,queries, aligns, targetIds, queryIds); break;
    case 24: FindPossibleOrientedAlignments<24>
               (target,queries, aligns, targetIds, queryIds); break;
    case 32: FindPossibleOrientedAlignments<32>
               (target,queries, aligns, targetIds, queryIds); break;
    case 48: FindPossibleOrientedAlignments<48>
               (target,queries, aligns, targetIds, queryIds); break;
    case 64: FindPossibleOrientedAlignments<64>
               (target,queries, aligns, targetIds, queryIds); break;
    case 96: FindPossibleOrientedAlignments<96>
               (target,queries, aligns, targetIds, queryIds); break;
    default:
      FatalErr("bad kmer_size: " << KMER_SIZE);
  }
}

