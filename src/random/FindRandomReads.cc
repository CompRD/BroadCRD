/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


#include "random/FindRandomReads.h"
#include "random/Random.h"

/// Generates positions randomly distributed on the specified
/// reference. Each generated position will be represented
/// as a pair (contig, offset_within_contig), so that, e.g.,
/// the base at the selected position is ref[contig][offset_within_contig].
/// Offset of the generated positions within their contigs will always be 
/// at least LEFT_MARGIN and will never exceed contig_size-RIGHT_MARGIN-1.
/// This method is safe in a sense that short contigs, on which no positions 
/// can be placed (i.e. contig_size < LEFT_MARGIN+RIGHT_MARGIN ) will be 
/// skipped and can never be returned; however, if there is no single contig
/// that is long enough to place a position onto, the vector \c positions
/// \em WILL have length 0 upon return, not the requested length \c n !
void FindRandomReads::Positions(int n, const vecbasevector & ref,
			vec<pair<int,int> > & positions,
			int LEFT_MARGIN ,
			int RIGHT_MARGIN ) {
    vec<int> sizes(ref.size());
    for (size_t i=0; i != ref.size(); ++i) {
      sizes[i] = ref[i].size();
    }
    Positions(n, sizes, positions, LEFT_MARGIN, RIGHT_MARGIN);
}

/// Same as the overloaded implementation 
/// (see Positions(int, const vecbasevector &, vec<pair<int,int> >, int, int)
/// for the details), but takes the vector of contig sizes as its argument
/// instead of the reference itself. Will return empty vector of positions
/// if all contigs are too short to place a position!
void FindRandomReads::Positions(int n, const vec<int> & sizes,
				vec<pair<int,int> > & positions,
				int LEFT_MARGIN,
				int RIGHT_MARGIN) {
  longlong totalsize=0;
  int excluded = RIGHT_MARGIN + LEFT_MARGIN;
  for (int i=0; i != sizes.isize(); ++i) {
    if ( excluded > sizes[i] ) continue; // contig i is too short!
    // add the number of allowed bases on contig i (those that
    // can be chosen):
    totalsize += ( sizes[i] - excluded ); 
  }
  for (int i=0; i != n; ++i) {
    longlong pos = longlong(drand48() * totalsize);
    for ( int j=0 ; j != sizes.isize() ; ++j ) {
      if ( excluded > sizes[j] ) continue; // skip contigs that are too short
      if ( pos < sizes[j] - excluded ) {
	positions.push_back( make_pair(j,pos + LEFT_MARGIN) );
        break;
      } else {
	pos -= ( sizes[j]- excluded );
      }
    }
  }
}

/// Generates \c n reads of length READ_SIZE randomly placed on 
/// the reference \c ref. If  \c reverseHalf is \c true (default),
/// each reads will be reverse complemented with 50% probability.
/// NOTE: if all the reference contigs are too short to accomodate a 
/// read, the returned vector \c reads will be empty! 
void FindRandomReads::Reads(int n, const vecbasevector & ref,
		    vecbasevector & reads, int READ_SIZE,
		    bool reverseHalf ) {
    vec<pair<int,int> > positions;
    Positions(n, ref, positions, 0, READ_SIZE);
    n = positions.size() ; // if returned positions vector is empty, it's now taken care of
    reads.Reserve( (n*READ_SIZE)/16 + n, n );
    basevector b;
    for (int i=0; i != n; ++i ) {
      b.SetToSubOf(ref[positions[i].first],positions[i].second,READ_SIZE);
      if (reverseHalf && (randomx() % 2)) b.ReverseComplement();
      reads.push_back(b);
    }
    
}   

