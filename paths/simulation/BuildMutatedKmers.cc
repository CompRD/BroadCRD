/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


// A class to track what part of the total 4^k space of kmers is
// actually occupied by a given set of basevectors.

#include "Basevector.h"
#include "kmers/TotalKmerSet.h"
#include "MainTools.h"
#include "ReadLocation.h"
#include "random/Random.h"

#include <set>

int main( int argc, char** argv ) 
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String( GENOME );
  CommandArgument_String( SPACE );
  CommandArgument_String( BAD_READS_PREFIX );
  CommandArgument_Int( COVERAGE );
  CommandArgument_Int( READ_LEN );
  CommandArgument_Double( ERR_RATE );
  CommandArgument_UnsignedInt_OrDefault( random_seed, 0 );
  EndCommandArguments;

  const int K = 20;
  TotalKmerSet<20> kmerSet;

  vecbasevector genome( GENOME );
    
  longlong genomeSize = 0;
  for ( size_t i = 0; i < genome.size(); ++i )
    genomeSize += genome[i].size();

  longlong numReads = ( genomeSize * (longlong) COVERAGE ) / (longlong)READ_LEN;

  PRINT( numReads );
  
  const int maxErrs = 4;
  double fracWithNumErrs[1+maxErrs];
  longlong countWithNumErrs[1+maxErrs];

  for ( int numErrs = 0; numErrs <= maxErrs; ++numErrs ) {
    double &prob = fracWithNumErrs[numErrs];
    prob = 1.0;
    for ( int i = numErrs; i < READ_LEN; ++i )
      prob *= ( 1.0 - ERR_RATE );
    for ( int i = 0; i < numErrs; ++i )
      prob *= ERR_RATE;
    for ( int i = 0; i < numErrs; ++i )
      prob *= (double)(READ_LEN-i) / (double)(i+1);

    countWithNumErrs[numErrs] = (longlong) round( prob * double( numReads ) );
    cout << "Reads with " << numErrs << " errors: " << countWithNumErrs[numErrs] << endl;
  }

  if ( ! IsRegularFile( SPACE ) ) {
    kmerSet.Build( genome );
    
    kmerSet.Write( SPACE );
  }
  else {
    kmerSet.Read( SPACE );
  }

  srandomx( random_seed );

  vecbasevector badReads;
  vec<read_location> badReadLocs;

  for ( int numErrs = 1; numErrs <= maxErrs; ++numErrs ) {
    longlong numReadsKept = 0;

    for ( int i = 0; i < countWithNumErrs[numErrs]; ++i ) {
      int genomeSeqId = randomx() % genome.size();
      int genomeSeqPos = big_random() % ( genome[genomeSeqId].size() - (READ_LEN-1) );
      
      basevector badRead;
      badRead.SetToSubOf( genome[genomeSeqId], genomeSeqPos, READ_LEN );

      int rc = ( randomx() % 2 );
      if ( rc == 1 )
        badRead.ReverseComplement();

      vector<bool> errPositions( READ_LEN, false );
      int errNum = 1;
      while ( errNum <= numErrs ) {
        int errPos = randomx() % READ_LEN;
        if ( ! errPositions[errPos] ) {
          int errAmt = 1 + randomx() % 3;
          badRead.Set( errPos, ( badRead[ errPos ] + errAmt ) % 4 );
          errPositions[errPos] = true;
          ++errNum;
        }
      }
      
      bool hasNonGenomicKmer = false;
      const int numKmers = READ_LEN - (K-1);
      for ( int kmerPos = 0; kmerPos < numKmers; ++kmerPos ) {
        basevector kmer;
        kmer.SetToSubOf( badRead, kmerPos, K );
        if ( ! kmerSet.Exists( kmer ) ) {
          hasNonGenomicKmer = true;
          break;
        }
      }

      if ( hasNonGenomicKmer )
        continue;

      badReadLocs.push_back( read_location( badReads.size(), READ_LEN, 
                                            genomeSeqId, genomeSeqPos,
                                            ( rc==1 ? ReverseOr : ForwardOr ),
                                            genome[genomeSeqId].size() ) );
      badReads.push_back_reserve( badRead );
      ++numReadsKept;
    }

    cout << "Reads kept with " << numErrs << " errors: " << numReadsKept << endl;
  }

  badReads.WriteAll( BAD_READS_PREFIX + ".fastb" );
  WriteLocs( BAD_READS_PREFIX + ".locs", badReadLocs );
}
  
