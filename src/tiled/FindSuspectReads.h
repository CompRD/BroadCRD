// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

// Given a set of read tuples between which a SNP exists,
// FindSnpCycles() finds 3-cycles (where read A has a SNP with read B,
// which has a SNP with read C, which has a SNP with read A, i.e., in
// terms of haplotype equality, A != B != C != A ) and 5-cycles (where
// read A has a SNP with read B, which is partnered with read B',
// which has a SNP with read C', which is partnered with read C, which
// has a SNP with read A, i.e., A != B == B' != C' == C != A).

// FindSuspectReads() returns in order the ids of the reads that
// participate in the most cycles.  The "most popular" read is found
// and noted, then all of the cycles in which it participates are
// removed, and the new "most popular" read is found and noted, and so
// on, until no cycles remain.

#include "ReadPairing.h"
#include "Vec.h"

#include <set>
#include <iostream>




struct snp_cycle
{
  snp_cycle() { }
    
  snp_cycle( set<int> &setReadsInCycle )
      : mSetReadsInCycle( setReadsInCycle )
  { }
    
  set<int> mSetReadsInCycle;

  friend bool operator==(const snp_cycle &lhs, const snp_cycle &rhs)
  {
      return ( lhs.mSetReadsInCycle == rhs.mSetReadsInCycle );
  }

  friend bool operator!=(const snp_cycle &lhs, const snp_cycle &rhs)
  {
      return !(lhs == rhs);
  }

  friend bool operator<(const snp_cycle &lhs, const snp_cycle &rhs)
  {
      pair<set<int>::const_iterator,set<int>::const_iterator> first_mismatch;
      first_mismatch = mismatch( lhs.mSetReadsInCycle.begin(), lhs.mSetReadsInCycle.end(),
                                 rhs.mSetReadsInCycle.begin() );

      if ( first_mismatch.first == lhs.mSetReadsInCycle.end() )
          if ( first_mismatch.second == rhs.mSetReadsInCycle.end() )
              return false;
          else
              return true;
      else
          if ( first_mismatch.second == rhs.mSetReadsInCycle.end() )
              return false;
          else
              return *(first_mismatch.first) < *(first_mismatch.second);
  }


  friend ostream& operator<<( ostream &ostrm, snp_cycle &sc )
  {
      copy( sc.mSetReadsInCycle.begin(), sc.mSetReadsInCycle.end(),
            ostream_iterator<int>( ostrm, " " ) );
      return ostrm;
  }

};


// if the contig is not consistent with two haplotypes, find the reads
// involved in the snafoo.

void FindSnpCycles( const vec< pair<int,int> > &snps, 
                    const vec<read_pairing> &pairs,
                    const vec<int> &pairs_index,
                    vec<snp_cycle> &vecCycles );

// in: a list of each set of reads that is inconsistent with the 2-hap
//     model for a single contig
// out: a list of read ids to remove from the contig in an attempt to get it to 
//     pass hap-consistency.
//
// idea: identify read(s) which are involved in the most inconsistent groups.
// note: in the case of a tie, simplistically chooses the read with the lowest id

vec<int> FindSuspectReads( const vec<snp_cycle> &snpCycles,
                           ostream *logp = 0 );
