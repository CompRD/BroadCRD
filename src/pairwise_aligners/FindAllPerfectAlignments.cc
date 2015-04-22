/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// This utility generates ALL the perfect alignments between sequences
// in some fastb file, or between sequences in two different fastb files.

#include "MainTools.h"

#include "pairwise_aligners/PerfectAligner.h"

// TODO: potentially dangerous truncation of index for all ints in PassInfo
typedef
struct {
  int firstBegin;
  int firstEnd;
  int secondBegin;
  int secondEnd;
  int partition;
  bool swap;
} PassInfo;

// Determine the number of passes needed for a MULTIPASS alignment,
// and all the requisite info.
void BuildPasses( const vecbasevector& sequences,
                  const int partition, 
                  vec<PassInfo>& passInfos );


int main( int argc, char **argv )
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String( FASTB1 );
  CommandArgument_String_OrDefault( FASTB2, "" );
  CommandArgument_String( OUT_ALIGNS );
  CommandArgument_UnsignedInt_OrDefault( PARTITION, 0 );
  CommandArgument_UnsignedInt_OrDefault( K, 24 );
  CommandArgument_Bool_OrDefault( INDEX, False );
  CommandArgument_Bool_OrDefault( FIND_SEMIPROPER, True );
  CommandArgument_Bool_OrDefault( FIND_IMPROPER, True );
  CommandArgument_Bool_OrDefault( MULTIPASS, False );
  CommandArgument_UnsignedInt_OrDefault( MAXKMERFREQ, 0 );
  EndCommandArguments;

  if ( ! FASTB2.empty() && PARTITION != 0 ) {
    cout << "Either specify FASTB2 or provide a non-zero PARTITION, but not both." << endl;
    TracebackThisProcess();
  }

  vec<alignment_plus> perfectAligns;

  PerfectAligner::Behavior desiredBehavior = PerfectAligner::findProperOnly;
  if ( FIND_SEMIPROPER )
    desiredBehavior = PerfectAligner::findSemiproper;
  if ( FIND_IMPROPER )
    desiredBehavior = PerfectAligner::findImproper;

  cout << Date() << ": Loading data." << endl;
  vecbasevector sequences( FASTB1 );
  cout << Date() << ": Done." << endl;
  
  int partition = PARTITION;
  if ( partition == 0 )
    partition = -1;

  if ( ! FASTB2.empty() && FASTB1 != FASTB2 )
  {
    partition = sequences.size();
    sequences.ReadRange( FASTB2, 0, MastervecFileObjectCount( FASTB2 ) );
  }
  
  if ( ! MULTIPASS ) {
    cout << Date() << ": Finding alignments." << endl;
    PerfectAligner aligner( K, desiredBehavior, &cout );
    aligner.SetMaxKmerFreq( MAXKMERFREQ );
    aligner.Align( sequences, perfectAligns, partition );
    cout << Date() << ": Done." << endl;
  }

  else {
    vec<PassInfo> passInfos;

    BuildPasses( sequences, partition, passInfos );

    vec<alignment_plus> passAligns;
    vecbasevector passSequences;
    passSequences.reserve( sequences.size() );

    cout << Date() << ": Aligning in " << passInfos.size() << " passes." << endl;

    for ( unsigned int pass = 0; pass < passInfos.size(); ++pass ) {
      cout << Date() << ": Pass " << pass+1 << endl;

      PassInfo& passInfo = passInfos[pass];

      cout << Date() << ":  Loading sequences." << endl;
      vec<Bool> seqsToCopy( sequences.size(), False );
      for ( int i = passInfo.firstBegin; i < passInfo.firstEnd; ++i )
        seqsToCopy[i] = True;
      for ( int i = passInfo.secondBegin; i < passInfo.secondEnd; ++i )
        seqsToCopy[i] = True;
      
      passSequences.clear();
      for ( size_t i = 0; i < sequences.size(); ++i )
        if ( seqsToCopy[i] )
          passSequences.push_back( sequences[i] );
        else
          passSequences.push_back( basevector() );

      cout << Date() << ":  Done." << endl;
      
      cout << Date() << ":  Finding alignments." << endl;
      passAligns.clear();

      PerfectAligner aligner( K, desiredBehavior, &cout );
      aligner.SetMaxKmerFreq( MAXKMERFREQ );
      aligner.Align( passSequences, passAligns, passInfo.partition );
      
      if ( ! passInfo.swap ) {
        copy( passAligns.begin(), passAligns.end(),
              back_inserter( perfectAligns ) );
      }
      else {
        for ( unsigned int i = 0; i < passAligns.size(); ++i ) {
          alignment_plus &thisAlign = passAligns[i];

          thisAlign.SafeSetId2( passInfo.partition + thisAlign.Id2() );

          perfectAligns.push_back( thisAlign );

          if ( thisAlign.Id1() != thisAlign.Id2() ) {
            int len1 = passSequences[ thisAlign.Id1() ].size();
            int len2 = passSequences[ thisAlign.Id2() ].size();
            thisAlign.Swap( len1, len2 );

            perfectAligns.push_back( thisAlign );
          }
        }
      }

      cout << Date() << ":  Done." << endl;
    }
  }

  cout << Date() << ": Sorting alignments." << endl;
  sort( perfectAligns.begin(), perfectAligns.end() );
  cout << Date() << ": Done." << endl;

  cout << Date() << ": Writing alignments." << endl;
  if ( INDEX )
    WriteAlignsWithIndex( OUT_ALIGNS, perfectAligns );
  else
  {
    ofstream out( OUT_ALIGNS.c_str() );
    out << perfectAligns;
  }
  cout << Date() << ": Done." << endl;

  return 0;
}


void BuildPasses( const vecbasevector& sequences, 
                  const int partition, 
                  vec<PassInfo>& passInfos ) {

  if ( partition == 0 ) {
    // If it's a self-compare, split the data in half by size; call
    // the first half "A" and the second half "B".  Do three passes: A
    // vs A, A vs B (adding swapped copies), and B vs B.

    passInfos.resize( 3 );

    size_t passSize = sequences.sumSizes()/2;

    size_t rawSizeSoFar = 0;
    size_t splitPoint = 0;
    for ( ; splitPoint < sequences.size(); ++splitPoint ) {
      rawSizeSoFar += sequences[splitPoint].size();
      if ( rawSizeSoFar >= passSize )
        break;
    }

    // AA
    passInfos[0].firstBegin = 0;
    passInfos[0].firstEnd = splitPoint;
    passInfos[0].secondBegin = 0;
    passInfos[0].secondEnd = splitPoint;
    passInfos[0].partition = 0;
    passInfos[0].swap = false;

    // AB
    passInfos[1].firstBegin = 0;
    passInfos[1].firstEnd = splitPoint;
    passInfos[1].secondBegin = splitPoint;
    passInfos[1].secondEnd = sequences.size();
    passInfos[1].partition = splitPoint;
    passInfos[1].swap = true;

    // BB
    passInfos[2].firstBegin = splitPoint;
    passInfos[2].firstEnd = sequences.size();
    passInfos[2].secondBegin = splitPoint;
    passInfos[2].secondEnd = sequences.size();
    passInfos[2].partition = 0;
    passInfos[2].swap = false;
  } 
  else {
    // If it's not a self-compare, find the sizes of the two
    // partitions.
    longlong firstPartTotal = 0;
    for ( int i = 0; i < partition; ++i ) 
      firstPartTotal += sequences[i].SizeOfDynamicData();
    
    longlong secondPartTotal = 0;
    for ( size_t i = partition; i < sequences.size(); ++i )
      secondPartTotal += sequences[i].SizeOfDynamicData();

    if ( firstPartTotal > secondPartTotal / 2 ) {
      if ( secondPartTotal > firstPartTotal / 2 ) {
        // Neither part is less than half the size of the other.
        passInfos.resize( 1 );

        // Compare all of the first part against all of the second part.
        passInfos[0].firstBegin = 0;
        passInfos[0].firstEnd = partition;
        passInfos[0].secondBegin = partition;
        passInfos[0].secondEnd = sequences.size();
        passInfos[0].partition = partition;
        passInfos[0].swap = false;
      }
      else {
        // The second part is less than half the size of the first.
        passInfos.resize( 2 );

        longlong rawSizeSoFar = 0;
        int splitPoint = 0;
        for ( ; splitPoint < partition; ++splitPoint ) {
          rawSizeSoFar += sequences[splitPoint].SizeOfDynamicData();
          if ( rawSizeSoFar >= firstPartTotal/2 )
            break;
        }

        // Compare the first half of the first part against all of the second part.
        passInfos[0].firstBegin = 0;
        passInfos[0].firstEnd = splitPoint;
        passInfos[0].secondBegin = partition;
        passInfos[0].secondEnd = sequences.size();
        passInfos[0].partition = partition;
        passInfos[0].swap = false;

        // Compare the second half of the first part against all of the second part.
        passInfos[1].firstBegin = splitPoint;
        passInfos[1].firstEnd = partition;
        passInfos[1].secondBegin = partition;
        passInfos[1].secondEnd = sequences.size();
        passInfos[1].partition = partition;
        passInfos[1].swap = false;
      }
    }
    else {
      // The first part is less than half the size of the second.
      passInfos.resize( 2 );
      
      longlong rawSizeSoFar = 0;
      size_t splitPoint = partition;
      for ( ; splitPoint < sequences.size(); ++splitPoint ) {
        rawSizeSoFar += sequences[splitPoint].SizeOfDynamicData();
        if ( rawSizeSoFar >= secondPartTotal/2 )
          break;
      }

      // Compare all of the first part against the first half of the second part.
      passInfos[0].firstBegin = 0;
      passInfos[0].firstEnd = partition;
      passInfos[0].secondBegin = partition;
      passInfos[0].secondEnd = splitPoint;
      passInfos[0].partition = partition;
      passInfos[0].swap = false;

      // Compare all of the first part against the second half of the second part.
      passInfos[1].firstBegin = 0;
      passInfos[1].firstEnd = partition;
      passInfos[1].secondBegin = splitPoint;
      passInfos[1].secondEnd = sequences.size();
      passInfos[1].partition = partition;
      passInfos[1].swap = false;
    }
  }
}
      
