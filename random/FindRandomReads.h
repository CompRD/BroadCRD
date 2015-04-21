/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

///Find random starting positions and/or sequences in a set of sequences.
/// \class FindRandomReads
///
/// Class for use in simulations where you need to pick random
/// reads or starting positions for reads/alignments from a genome.

#include "Basevector.h"
#include "CoreTools.h"
#include "random/Random.h"
#include "math/Functions.h"

class FindRandomReads {
public:
  /// Same as the overloaded implementation 
  /// (see Positions(int, const vecbasevector &, vec<pair<int,int> >, int, int)
  /// for the details), but takes the vector of contig sizes as its argument
  /// instead of the reference itself. Will return empty vector of positions
  /// if all contigs are too short to place a position!
  static void Positions(int n, const vec<int> & sizes,
			vec<pair<int,int> > & positions,
			int LEFT_MARGIN = 0,
			int RIGHT_MARGIN = 0);

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
  static void Positions(int n, const vecbasevector & ref,
			vec<pair<int,int> > & positions,
			int LEFT_MARGIN = 0,
			int RIGHT_MARGIN = 0);

  /// Generates \c n reads of length READ_SIZE randomly placed on 
  /// the reference \c ref. If  \c reverseHalf is \c true (default),
  /// each reads will be reverse complemented with 50% probability.
  /// NOTE: if all the reference contigs are too short to accomodate a 
  /// read, the returned vector \c reads will be empty! 
  static void Reads(int n, const vecbasevector & ref,
		    vecbasevector & reads, int READ_SIZE,
		    bool reverseHalf = True);

  // This version of Reads picks reads with a given first base distribution.  You 
  // should have Sum(first_base_distribution) = 1.0.

  static void Reads( int n, const vec<double>& first_base_distribution, 
          const vecbasevector & ref, vecbasevector & reads, int READ_SIZE,
          bool reverseHalf = True ) 
    {    
      // Reserve space.

      reads.Reserve( (n*READ_SIZE)/16 + n, n );

      // Get four base counts that add exactly to n.

      vec<int> counts(4, 0);
      vec<double> target(4);
      for ( int i = 0; i < 4; i++ )
	target[i] = int( round( first_base_distribution[i] * double(n) ) );
      for ( int i = 0; i < n; i++ )  {    
	double M = Max(target);
	int j = Position( target, M );
               ++counts[j];
               --target[j];    
      }

      // Get the reads.

      while( Sum(counts) > 0 ) {    
	vec<pair<int,int> > positions;
	Positions( n, ref, positions, 0, READ_SIZE );
	basevector b;
	for (int i=0; i != n; ++i ) {    
	  b.SetToSubOf(
		     ref[positions[i].first], positions[i].second, READ_SIZE );
	  if (reverseHalf && (randomx() % 2)) b.ReverseComplement();
	  int firstbase = b[0];
	  if ( counts[firstbase] == 0 ) continue;
	  --counts[firstbase];
	  reads.push_back(b);    
	}    
      }    
    }

};
