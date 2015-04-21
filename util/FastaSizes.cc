/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/// FastaSizes: print out the sizes (in bases) of the records in a fasta file.
/// \file FastaSizes.cc
/// Reads from standard input.  
/// Writes a single line to stdout, with sizes separated by blanks.

#include "MainTools.h"
  
inline void readLine() {
  static int c;
  while ( (c = getchar()) != '\n' && c != EOF) {};
}

int main()
{
  int c = getchar();
  if (c != '>') {
    printf("Not in fasta format!");
    return 1;
  }
  readLine();

  int count = 0;  
  while ((c = getchar()) != EOF) {
    if (c == '>') {
      readLine();
      printf("%d ", count);
      flush(cout);
      count = 0;
    }
    else if (isalpha(c)) ++count;
  }
  printf("%d\n", count);
  return 0;
}
  
