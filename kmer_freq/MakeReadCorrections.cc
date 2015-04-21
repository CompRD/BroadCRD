/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include <map>
#include "kmer_freq/MakeReadCorrections.h"
#include "kmer_freq/WriteKmerFrequencies.h"
#include "math/Functions.h"
#include "paths/BaseErrorProb.h"


bool findNextSeed( vec<bool> strongPositions, int currentRightEdge, 
                   int& newLeftEdge, int& newRightEdge )
{
  newLeftEdge = newRightEdge = currentRightEdge; // initialization
  for ( int pos = currentRightEdge + 1; pos < strongPositions.isize(); ++pos ) {
    if ( strongPositions[ pos ] ) {
      newLeftEdge  = pos;
      newRightEdge = pos;
      break;
    }
  }
  for ( int pos = newLeftEdge + 1; pos < strongPositions.isize(); ++pos ) {
    if( ! strongPositions[ pos ] ) { 
      newRightEdge = pos - 1; 
      break;
    }
  }
  if ( newLeftEdge > currentRightEdge )
    return true;
  else
    return false;
}  


/**
   \file

   Actual implementation of ProvisionalEdits.cc .

   \copydoc ProvisionalEdits.cc
   \ingroup grp_edits
*/

/**
   \copydoc ProvisionalEdits.cc

   This routine corrects the reads according to the procedure described above.
   
   \param[in,out] reads              the reads; any edits (fixes) will be applied in-place to the reads in this array.
   \param[in]     i_reads_todo       the indices of the reads to correct
   \param[in]     quals              the \link grp_quals quality scores\endlink of the reads
   \param[in]     tablePtr           the strong KmerShortMap
   \param[out]    read_is_strong     whether the read is determined to be strong or not
   \param[in]     verbose            whether to print more debug info
   \param[in]     probProfile        if quality scores not given in \p quals, then the probability
                                     of error at each read \e position
   \param[in]     max_errors         consider mutations with up to this many errors
*/
void MakeReadCorrections(vecbasevector& reads, 
                         const vec<int>& i_reads_todo,
                         const vecqualvector& quals,
                         const KmerShortMap& KmerMap,
                         vec<Bool>& read_is_strong,
                         const bool verbose,
                         const int max_errors,
                         const bool trim_reads,
                         const bool thorough,
			 const bool keep_partial)
{

  
  // Make changes.
  cout << Date() << " starting making changes" << endl;
  longlong numchanged             = 0;
  longlong numweak                = 0;
  longlong numStrongInitialSeeds  = 0;
  longlong numStrongModifiedSeeds = 0;
  longlong numExtensionsAttempted = 0;
  longlong sumReadLens            = 0;
  longlong sumCrtdReadLens        = 0;
  longlong singlemut              = 0;
  longlong multiplemut            = 0;
  longlong numberOfTrials         = 0;
  longlong numStrongFromModified  = 0;

  // Find minimum and maximum value of K
  int Ksize    = KmerMap.GetKmerSize();
  int seedSize = Ksize; /// the size of a seed used for extension of strong k-mers

  // If error profile not passed then create default flat profile for largest read size
  BaseErrorProbProfile probProfile = BaseErrorProbProfile(seedSize);

  BaseErrorTable  seedErrorTable = probProfile.getErrorTable(seedSize, 1, 0, false);
  size_t n_todo = i_reads_todo.size();
  double clock = -1.0;
  for (size_t i_todo = 0; i_todo < n_todo; i_todo++) {
    int i_read = i_reads_todo[i_todo];
    read_is_strong[i_read] = False;  // reads are weak by default
    



    if (i_todo % 20000 == 0) {
      double time_used = WallClockTime( ) - clock;
      cout << "  " << Date() 
           << " i_todo = " << i_todo << " of " << n_todo
           << ", num weak  = "   << numweak
           << ", num changed = " << numchanged
           << ", time used = "   << setprecision(3) << time_used << " seconds" 
           << ", reads left " << (n_todo - i_todo - 1) 
           << endl;
      clock = WallClockTime( );    
    }
    
    if ( verbose ) 
      cout << "\n\n--- NEW READ i_todo = " << i_todo << " ---" << endl;




    basevector& theRead = reads[i_read];
    
    if ( KmerMap.IsStrong( theRead ) ) {  
      read_is_strong[i_read] = True; // read is strong
      if ( verbose ) 
        cout << "READ i_read = " << i_read << " is strong originally.\n";
      continue;
    }
    ++numweak; // read is not strong

    int readSize = theRead.isize();

    /// Cannot correct reads shorter than the seed size in tables
    if (readSize < seedSize) {
      if ( verbose ) cout << "this read is too short, skipping" << endl;
      continue;
    }

    basevector mutatedRead = theRead;

    /// Marking strong positions in the original read and looking for a starting position 
    /// for "strong" k-mer extension (choosing a leftmost one, presumably a higher quality part of a read)
    bool seedFound  = false;
    int  startPos   = 0;
    vec<int> kmerValues( theRead.isize(), 0 );
    vec<bool> strongPositions( theRead.isize(), false );
    KmerMap.GetValues( theRead, kmerValues );
    for ( int pos = 0; pos < theRead.isize() - Ksize +1; ++pos )
      if ( kmerValues[pos] > 0 ) {
        if ( ! seedFound ) {
          seedFound = true;
          startPos  = pos;
          ++numStrongInitialSeeds;
        }
        for ( int pos2 = pos; pos2 < pos + Ksize; ++pos2 )
          strongPositions[ pos2 ] = true;
      }

    if ( verbose ) {
      cout << "initial strong positions:" << endl;
      strongPositions.Println( cout );
    }
    bool origSeedFound = seedFound;
 
    /// ==================================================================================
    /// This code section finds strong seeds in the read by mutating the bases
    if ( !seedFound || thorough ) {
      if ( verbose )
        cout << "uncovering strong positions" << endl;
      // To hold record of best sucessful change
      vec<int> mutpos_hit;
      vec<int> mutbase_hit;
      // Current change information
      vec<int> mutpos;
      vec<int> mutbase;
      
      int seedMutTableSize = seedErrorTable.isize();
      int mutationno;
     
      bool localSeedFound = false;
      for ( int iStartPos = 0; iStartPos < theRead.isize() - seedSize + 1; ++iStartPos ) {
        if ( strongPositions[ iStartPos ] ) continue;
        if ( localSeedFound && !thorough ) break;
       
        int mutcountLast      = 0;
        int numStrongFound    = 0;
        for( mutationno = 0; mutationno < seedMutTableSize; ++mutationno ) {
          // Get potential change information
          if ( numStrongFound > 1 ) 
            break; /// to many alternatives
          mutpos       = seedErrorTable[ mutationno ].base_pos;
          int mutcount = mutpos.isize();
          if ( mutcountLast < mutcount ) {
            mutcountLast = mutcount;
            if ( numStrongFound == 1 ) 
              break; // was strong with fewer mutations, stop
          }
          mutbase.resize(mutcount);
          // Try out all possible base/position combinations
          int nocombs = 2 << mutcount * 2 - 1;
          for (int comb = 0; comb < nocombs; ++comb) {
            bool skip    = false;
            int basesum  = comb;
            for (int posno = 0; posno < mutcount; ++posno) {
              int base = basesum % 4;
              int pos  = mutpos[posno] + iStartPos;
              if ( strongPositions[ pos ] ) {
                skip = true;
                break;
              }
              else {
                if (theRead[pos] == base) {
                  skip = true;
                  break;
                } 
                else {
                  mutbase[posno] = base;
                  mutatedRead.Set( pos, base );
                  basesum /= 4;
                }
              }
            } 
            if (skip) 
              continue; // Current base/position combination same as original - skip. 
	    
            numberOfTrials++;
            // Measure effect of base changes.
            if ( KmerMap.IsStrong( mutatedRead, iStartPos, iStartPos + seedSize -1 )  ) {
              ++numStrongFound;
              if ( numStrongFound > 1 ) 
                break;
              // Found a sucessful change - record it and try another combination.
              mutpos_hit        = mutpos;
              mutbase_hit       = mutbase;
            }
          }
          // reset to the original base before moving to the next position
          for (int posno = 0; posno < mutcount; ++posno)
            mutatedRead.Set( mutpos[posno] + iStartPos, theRead[ mutpos[posno] + iStartPos ] );
        }
        // Modify read if a single sucessful change was found
        if ( numStrongFound == 1 ) {
          localSeedFound = true;
          numStrongModifiedSeeds++;
          for (int posno = 0; posno < mutpos_hit.isize(); ++posno)
            mutatedRead.Set( mutpos_hit[posno] + iStartPos, mutbase_hit[posno] );
          for ( int ipos = iStartPos; ipos < iStartPos + seedSize; ++ipos )
            strongPositions[ ipos ] = true;
          iStartPos += seedSize -1;
        }
      }
      
      if ( verbose ) {
        cout << "uncovered strong positions :" << endl;
        strongPositions.Println( cout );
      }
      
      seedFound          = false;
      startPos           = 0;
      for ( int pos = 0; pos < strongPositions.isize(); ++pos ) {
        if ( strongPositions[ pos ] ) {
          seedFound = true;
          startPos  = pos;
          break;
        }
      }

    }
    
    if ( ! seedFound ) {
      if ( verbose ) cout << "did not find seeds in this read! Going to the next one\n" << endl;
      continue; // no point in going further since the seed was not found
    }
    
    
    ++numExtensionsAttempted;
    bool BadRead    = false;  // the read will be bad if it contains inconsistent strong kmers
                              // or too many mutations etc.
    // maximally extend the seeds 
    if ( verbose ) cout << "extending strong seeds" << endl;
    int lastRightPos, newLeftPos, newRightPos; 
    vec<int> matchValues( theRead.isize(), 0 );
    vec<int> alignmentValues( theRead.isize(), 0 );
    int matchScore     =  1;
    int mutationScore  = -1;
    int undefScore     = -2 * matchScore * theRead.isize(); // matches cannot offset this value  
    newLeftPos = newRightPos = lastRightPos = -1;
    while ( findNextSeed( strongPositions, lastRightPos, newLeftPos, newRightPos ) ) {
      if ( ! KmerMap.IsStrong( mutatedRead, newLeftPos, newRightPos ) ) { 
        BadRead = true;
        if ( verbose ) 
          cout << "inconsistent read" << endl;
        break;
      }
      for ( int k = newLeftPos; k <= newRightPos; ++k ) {
        if ( mutatedRead[ k ] == theRead[ k ] ) 
          matchValues[ k ] = matchScore;
        else 
          matchValues[ k ] = mutationScore;
      }

      while ( newLeftPos > 0 ) {
        --newLeftPos;
        if ( KmerMap.IsStrong( mutatedRead, newLeftPos, newLeftPos + Ksize -1 ) ) {
          strongPositions[ newLeftPos ]  = true;
          matchValues[ newLeftPos ]      = matchScore;
          continue;
        }
        else { 
          if ( strongPositions[ newLeftPos ] ) { // this position was determined as strong previously
            BadRead = true;
            ++newLeftPos;
            if ( verbose ) 
              cout << "inconsistent read" << endl;
            break;
          }
	  
          int nStrongLeftMut          = 0;
          int goodBase                = theRead[ newLeftPos ];
          for ( int base = 0; base < 4; ++base ) {
            if ( base == theRead[ newLeftPos ] ) 
              continue;
            mutatedRead.Set( newLeftPos, base );
            bool strongLeftMut = KmerMap.IsStrong( mutatedRead, newLeftPos, newLeftPos + Ksize -1 ) ? true : false;
            if ( strongLeftMut ) { 
              nStrongLeftMut++; 
              goodBase = base; 
            }
          }
          if ( nStrongLeftMut == 0 || nStrongLeftMut > 1 ) {
            /// there are none or too many strong mutations
            matchValues[ newLeftPos ]     = undefScore;
            break;
          }
          else {
            matchValues[ newLeftPos ]     =  mutationScore;
            strongPositions[ newLeftPos ] =  true;
            mutatedRead.Set( newLeftPos, goodBase );
          }
        }
      }
      
      if ( BadRead ) break;
      
   
      while ( newRightPos < readSize - 1 ) {
        ++newRightPos;
        if ( KmerMap.IsStrong( mutatedRead, newRightPos - Ksize +1, newRightPos ) ) {
          strongPositions[ newRightPos ]  = true;
          matchValues[ newRightPos ]      = matchScore;
          continue;
        }
        else { 
          if ( strongPositions[ newRightPos ] ) { // this position was deterined to be strong previously
            BadRead = true;
            --newRightPos;
            if ( verbose ) 
              cout << "inconsistent read" << endl;
            break;
          }
	  
          int nStrongRightMut         = 0;
          int goodBase                = theRead[ newRightPos ];
          for ( int base = 0; base < 4; ++base ) {
            if ( base == theRead[ newRightPos ] ) 
              continue;
            mutatedRead.Set( newRightPos, base );
            if ( KmerMap.IsStrong( mutatedRead, newRightPos - Ksize +1, newRightPos ) ) {
              nStrongRightMut++; goodBase = base; 
            }
          }
          if ( nStrongRightMut == 0 || nStrongRightMut > 1 ) {
            /// there are none or too many strong mutations
            matchValues[ newRightPos ] = undefScore;
            break;
          } 
          else {
            matchValues[ newRightPos ]     =  mutationScore;
            strongPositions[ newRightPos ] =  true;
            mutatedRead.Set( newRightPos, goodBase );
          }
        }
      }
      lastRightPos = newRightPos;
      if ( BadRead ) 
        break;
    }

    if ( BadRead ) 
      continue;

    if ( verbose ) {
      cout << "final strong positions:" << endl;
      strongPositions.Println( cout );
    }


    /// trim the read if requested  
    int finalLeftPos   = 0;
    int finalRightPos  = readSize -1;
    int countMutations = 0;
    if ( trim_reads ) {
      int aliValue      =  0;
      int bestAliValue  =  0;
      int bestPos       = -1;
      for ( int k = 0; k < matchValues.isize(); ++k ) {
        aliValue += matchValues[ k ];
        if ( aliValue <= 0 ) {
          alignmentValues[ k ] = 0;
        }
        else {
          alignmentValues[ k ] = aliValue;
          if ( aliValue > bestAliValue ) {
            bestAliValue  = aliValue;
            bestPos       = k;
          }
        }
      }
      // trace back
      if ( bestAliValue > 0 ) {
        finalRightPos = bestPos;
        for ( int k = finalRightPos -1; k >= 0; --k ) {
          if ( alignmentValues[ k ] <= 0 ) {
            finalLeftPos = k + 1;
            break;
          }
        }
        if ( KmerMap.IsStrong( mutatedRead, finalLeftPos, finalRightPos ) ) {
          theRead.SetToSubOf( mutatedRead, finalLeftPos, finalRightPos - finalLeftPos +1 );
          if ( ! origSeedFound ) 
            ++numStrongFromModified;
        }
        else {
          BadRead = true;
        }
      } 
      else {      
        BadRead = true; 
      }
    } 
    else if ( keep_partial ){
      theRead = mutatedRead;
      if ( ! origSeedFound && KmerMap.IsStrong( mutatedRead ) ) 
	++numStrongFromModified;
    }
    else{
      // check if the read is consistent and has not too many errors
      if ( !KmerMap.IsStrong( mutatedRead ) ) {
        BadRead = true;
      }
      else {
        for ( int k = 0; k < matchValues.isize(); ++k )
          if ( matchValues[ k ] != matchScore )
            ++countMutations;
        if ( countMutations > max_errors ) {
          BadRead = true;
        }
        else { 
          theRead = mutatedRead;
          if ( ! origSeedFound ) 
            ++numStrongFromModified;
        }
      }
    }

    if ( BadRead ) 
      continue;


    read_is_strong[i_read] = True;     // if we got here then the read is strong

    if ( verbose ) 
      cout << "READ i_read = " << i_read << " is strong\n";
    

    ++numchanged;
    if ( countMutations == 1 && finalLeftPos == 0 && finalRightPos == readSize -1 )
      ++singlemut;
    else
      ++multiplemut;
    
    sumReadLens      += readSize;
    sumCrtdReadLens  += theRead.isize();

    if ( verbose ) 
      cout << "read finished ---------\n\n";
  }
  PRINT2( numweak, numchanged );
  PRINT2( singlemut, multiplemut );
  cout << "numberOfTrials         = " << numberOfTrials << endl;
  cout << "numStrongInitialSeeds  = " << numStrongInitialSeeds << endl;
  cout << "numStrongModifiedSeeds = " << numStrongModifiedSeeds << endl;
  cout << "numExtensionsAttempted = " << numExtensionsAttempted << endl;
  cout << "numStrongFromModified  = " << numStrongFromModified << endl;
  cout << "sumReadLens            = " << sumReadLens << endl;
  cout << "sumCrtdReadLens        = " << sumCrtdReadLens << endl;
  double percNotTrimmed = 100.0 * sumCrtdReadLens / sumReadLens;
  cout << "percent read bases remaining after trimming = " << percNotTrimmed << endl;
}


template <class KSHAPE>
void MakeReadCorrections( vecbasevector& reads,  
                          const vec<int>& i_reads_todo,
                          const vecqualvector& quals,
                          vec<Bool>& read_is_strong,
                          const String& filename,
                          const bool verbose,
                          const int max_errors,  
                          const bool trim_reads,
                          const bool thorough,
			  const bool keep_partial) 
{
  KmerShortMap KmerMap( KSHAPE::getId(), filename );
  
  MakeReadCorrections( reads, i_reads_todo, quals, KmerMap, read_is_strong, 
                       verbose, max_errors );
}

#define INSTANTIATE(KSHAPE, dummy)                                      \
  template void MakeReadCorrections<KSHAPE>( vecbasevector&, const vec<int>&, \
                                             const vecqualvector&, vec<Bool>&, \
                                             const String&, const bool, const int, \
                                             const bool, const bool, const bool)
FOR_ALL_KSHAPES(INSTANTIATE,unused);
