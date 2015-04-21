///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/**
   Compare two fastb files and report any differences.
*/

const char *DOC =
"Compare two fastb files and report any differences. "
"Each basevector in the first file is compared with the basevector at the "
"corresponding index in the second file. If SORT=True then both vecbasevectors "
"are first sorted before being compared. In this way you can determine if "
"two fastb files contain the same set of basevectors if their orders differ.";


#include "MainTools.h"
#include "Basevector.h"
#include "math/Functions.h"

int main( int argc, char *argv[] ) 
{
  RunTime( );

  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String_Doc(IN1, "First fastb file.");
  CommandArgument_String_Doc(IN2, "Second fastb file.");
  CommandArgument_Bool_OrDefault_Doc(CANONICALIZE, False,
    "Canonicalize basevectors before comparison (and sorting).");
  CommandArgument_Bool_OrDefault_Doc(SORT, False,
    "Sort basevectors before comparison.");
  CommandArgument_Bool_OrDefault_Doc(BRIEF, False,
    "Only report if files differ.");
  CommandArgument_Bool_OrDefault_Doc(COUNT, False,
    "Only report count of differences.");
  CommandArgument_Bool_OrDefault_Doc(FAIL, False,
    "Return 1 if files differ.");
  EndCommandArguments;

  vecbasevector fastb1(IN1);
  vecbasevector fastb2(IN2);

  int n1 = fastb1.size();
  int n2 = fastb2.size();

  int i = 0;
  int diffCount = 0;

  if (CANONICALIZE) {
    for ( size_t i = 0; i < fastb1.size(); i++)
      fastb1[i].Canonicalize();
    for ( size_t i = 0; i < fastb2.size(); i++)
      fastb2[i].Canonicalize();
  }

  if (SORT) {
    fastb1.Sort();
    fastb2.Sort();
  }
  
  while (i < n1 && i < n2) {
    if (fastb1[i] != fastb2[i]) {
      if (!COUNT && !BRIEF)
	cout << "diff at index " << i << "\n";
        if ( fastb1[i].size( ) != fastb2[i].size( ) )
             cout << "different sizes\n";
        else
        {    cout << "same sizes\n";
             int diffs = 0;
             for ( size_t j = 0; j < fastb1[i].size( ); j++ )
             {    if ( fastb1[i][j] != fastb2[i][j] )
                  {    if ( diffs == 0 ) 
                            cout << "first base where they differ is " << j << "\n";
                       diffs++;    }    }    
             cout << "total substitutions = " << diffs << "\n";    }
      diffCount++;
      if (BRIEF)
	break;
    }
    ++i;
  }

  diffCount += Abs(n1 - n2);

  if (BRIEF && diffCount != 0) {
    cout << "Files differ\n";
  } else if (COUNT) {
    cout << diffCount << "\n";
  }
  
  if ( FAIL && diffCount > 0 ) return 1;
  return 0;
}
