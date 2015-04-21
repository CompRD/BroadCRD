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

int main( int argc, char** argv ) 
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String( FASTB );
  CommandArgument_String( SPACE );
  EndCommandArguments;

  const int K = 20;
  TotalKmerSet<20> kmerSpace;
  
  vecbasevector seqs( FASTB );
    
  kmerSpace.Build( seqs );
  kmerSpace.Write( SPACE );
}
  
