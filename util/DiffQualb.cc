/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   Compare two qualb files and report any differences.
*/

const char *DOC =
"Compare two qualb files and report any differences. "
"Each basevector in the first file is compared with the basevector at the "
"corresponding index in the second file. If SORT=True then both vecbasevectors "
"are first sorted before being compared. In this way you can determine if "
"two qualb files contain the same set of basevectors if their orders differ.";


#include "MainTools.h"
#include "Qualvector.h"
#include "math/Functions.h"

int main( int argc, char *argv[] ) 
{
  RunTime( );

  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String_Doc(IN1, "First qualb file.");
  CommandArgument_String_Doc(IN2, "Second qualb file.");
  CommandArgument_Bool_OrDefault_Doc(SORT, False,
    "Sort qualvectors before comparison.");
  CommandArgument_Bool_OrDefault_Doc(BRIEF, False,
    "Only report if files differ.");
  CommandArgument_Bool_OrDefault_Doc(COUNT, False,
    "Only report count of differences.");
  EndCommandArguments;

  vecqualvector qualb1(IN1);
  vecqualvector qualb2(IN2);

  int n1 = qualb1.size();
  int n2 = qualb2.size();

  PRINT2(n1,n2);

  int i = 0;

  int diffCount = 0;

  if (SORT) {
    qualb1.Sort();
    qualb2.Sort();
  }
  
  while (i < n1 && i < n2) {
    if (qualb1[i] != qualb2[i]) {
      if (!COUNT && !BRIEF)
	cout << "diff at index " << i << "\n";
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
  
  return 0;
}
