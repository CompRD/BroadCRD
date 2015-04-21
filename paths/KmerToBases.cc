// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

// KmerToBases.cc: interactive tool to query a KmerBaseBroker.

#include "MainTools.h"
#include "paths/KmerBaseBroker.h"

int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_UnsignedInt_OrDefault(K, 48);
  EndCommandArguments;

  String run_dir = PRE + "/" + DATA + "/" + RUN;

  KmerBaseBroker kbb( run_dir, K);  // true = verbose

  longlong k1, k2;
  basevector kmer_bases;

  while( cin ) {
    cout << "kmers: " << flush;
    cin >> k1 >> k2;
    
    if( !cin ) {
      cout << endl;
      break;
    }
    
    int d = kbb.MinOffset(k1,k2);
    cout << "Minimum offset " << d << endl;

    kbb.Bases(k1).Print(cout);
    while(d--) cout << ' ';
    kbb.Bases(k2).Print(cout);
    cout << endl;
  }

  EXIT_MAIN_NORMALLY;
}

