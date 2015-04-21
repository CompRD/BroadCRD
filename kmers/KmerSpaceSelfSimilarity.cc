/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


#include "Basevector.h"
#include "kmers/TotalKmerSet.h"
#include "MainTools.h"
#include "ReadLocation.h"
#include "random/Random.h"

#include <set>

template <int K>
void DoIt( const String& GENOME, const int max_muts, const int trials_per_mut ) {
  TotalKmerSet<K> kmerSpace;

  vecbasevector genome( GENOME );
  kmerSpace.Build( genome );

  for ( int numMuts = 1; numMuts <= max_muts; ++numMuts ) {
    int mutHits = 0;
    for ( int i = 0; i < trials_per_mut; ++i ) {
      int genomeSeqId = randomx() % genome.size();
      int genomeSeqPos = big_random() % ( genome[genomeSeqId].size() - (K-1) );
      
      basevector mutatedKmer;
      mutatedKmer.SetToSubOf( genome[genomeSeqId], genomeSeqPos, K );

      vector<bool> errPositions( K, false );
      int errNum = 1;
      while ( errNum <= numMuts ) {
        int errPos = randomx() % K;
        if ( ! errPositions[errPos] ) {
          int errAmt = 1 + randomx() % 3;
          mutatedKmer.Set( errPos, ( mutatedKmer[ errPos ] + errAmt ) % 4 );
          errPositions[errPos] = true;
          ++errNum;
        }
      }
      
      if ( kmerSpace.Exists( mutatedKmer ) )
        ++mutHits;
    }

    cout << numMuts << ": " << (double) mutHits / (double) trials_per_mut * 100.0 << "%" << endl;
  }
}  

int main( int argc, char** argv ) 
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String( GENOME );
  CommandArgument_Int( max_muts );
  CommandArgument_Int( trials_per_mut );
  CommandArgument_Int( K );
  CommandArgument_UnsignedInt_OrDefault( random_seed, 0 );
  EndCommandArguments;

  srandomx( random_seed );
  
  switch ( K ) {
    case 12:
      DoIt<12>( GENOME, max_muts, trials_per_mut );
      break;
    case 16:
      DoIt<16>( GENOME, max_muts, trials_per_mut );
      break;
    case 20:
      DoIt<20>( GENOME, max_muts, trials_per_mut );
      break;
    case 24:
      DoIt<24>( GENOME, max_muts, trials_per_mut );
      break;
    case 32:
      DoIt<32>( GENOME, max_muts, trials_per_mut );
      break;
  }
}
  
