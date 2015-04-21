// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
//
// Compare outputs of NamesToProject.cc from two runs, one with 'origDir' specified, one not.
// If 'origDir' is specified, read data from 'origDir'/{reads_orig.fastn, reads_orig.qualb,
// reads.traceinfo.tmp} instead of parsing the original data files.  Runs significantly faster 
// ( ~20 mins on imouse vs a few hours ), but also uses ~10Gb (on imouse)
//
// This test compares the sequence and quality scores for each read, exiting on single base
// differences or more broad differences if they exist (i.e. different number of reads 
// for each run)
//

#include "MainTools.h"
#include "Basevector.h"
#include "CompressedSequence.h"
#include "FastaFileset.h"
#include "Qualvector.h"

int main( int argc, char *argv[] )
{
  RunTime();

  BeginCommandArguments;
     CommandArgument_String(sequence_control);
     CommandArgument_String(sequence_test);
     CommandArgument_String(qualscores_control);
     CommandArgument_String(qualscores_test);
  EndCommandArguments;

  // first check the sequences

  cout << Date() << " Testing sequences." << endl;
  
  LastWordParser name_parser;
  FastaPairedFileset controlFileset( vec<String>(1,sequence_control),
                                     vec<String>(1,qualscores_control),
                                     &name_parser );

  FirstWordParser name_parser2;
  FastaPairedFileset testFileset( vec<String>(1,sequence_test),
                                  vec<String>(1,qualscores_test),
                                  &name_parser2 );

  vec<String> controlUnmatchedSequences, controlUnmatchedQuals;
  vec<String> testUnmatchedSequences, testUnmatchedQuals;

  controlFileset.GetUnmatchedSequenceNames( controlUnmatchedSequences );
  controlFileset.GetUnmatchedQualityNames( controlUnmatchedQuals );

  testFileset.GetUnmatchedSequenceNames( testUnmatchedSequences );
  testFileset.GetUnmatchedQualityNames( testUnmatchedQuals );

  if ( controlUnmatchedSequences.size() != testUnmatchedSequences.size() )
  {
    FatalErr( "Test set has " << testUnmatchedSequences.size() << " sequences without"
              << " quality scores, while control set has " 
              << controlUnmatchedSequences.size() << "." );
  }

  if ( controlUnmatchedQuals.size() != testUnmatchedQuals.size() )
  {
    FatalErr( "Test set has " << testUnmatchedQuals.size() << " quality scores without"
              << " sequences, while control set has " 
              << controlUnmatchedQuals.size() << "." );
  }

  String controlName, testName;
  CompressedSequence controlSeq, testSeq;
  qualvector controlQuals, testQuals;

  controlFileset.GetNext( controlName, controlSeq, controlQuals );
  testFileset.GetNext( testName, testSeq, testQuals );

  while ( 1 )
  {
    if ( controlName < testName )
    {
      cout << "Test set is missing sequence named: " << controlName << endl;
      if ( controlFileset.GetNext( controlName, controlSeq, controlQuals ) )
        continue;
      else
        break;
    }

    if ( testName < controlName )
    {
      cout << "Test set has extra sequence named: " << controlName << endl;
      if ( testFileset.GetNext( testName, testSeq, testQuals ) )
        continue;
      else
        break;
    }

    vec<char> controlChrs = controlSeq.asVecChar();
    vec<char> testChrs  = testSeq.asVecChar();

    for ( int j=0; j<(int) controlChrs.size(); ++j )
      if ( controlChrs[j] != testChrs[j] )
	FatalErr("Different bases at " << j <<" in " << controlName << ": "<< controlChrs[j]
		<< " " << testChrs[j] );

    for ( int j=0; j<(int) controlQuals.size(); ++j )
      if ( controlQuals[j] != testQuals[j] )
	FatalErr("Different quals at " << j <<" in "<< controlName << ": "<< controlQuals[j]
		<< " " << testQuals[j] );
  }

  cout << "\nVerified!\n" << endl;
}
