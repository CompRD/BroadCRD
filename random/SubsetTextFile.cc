/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "random/Shuffle.h"
#include "FastaFileset.h"

/** Pick LINES lines randomly from IN without replacement NUM_FILES times.
And put these sets of lines into NUM_FILES separate files.

\file SubsetTextFile.cc

*/

int main( int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String(IN); \
  CommandArgument_Int_OrDefault(LINES,80); 
  CommandArgument_Int_OrDefault(SEED,0); 
    CommandArgument_Int_OrDefault(NUM_FILES,100); 
  EndCommandArguments;


  Ifstream(is, IN);
  vec<String> file;
  String l;
  while (true) {
    getline(is, l);
    if (!is) break;
    file.push_back(l);
  }

  vec<int> shuffled;
  for (int i=0; i != NUM_FILES; ++i) {
    Shuffle(file.size(), shuffled, SEED+17*i);
    Ofstream(os, IN+"." + ToString(LINES)+"."+ToString(i) );
    for (int j=0; j != LINES; ++j) {
      os << file[shuffled[j]] << endl;
    }
  }

  return 0;
}

