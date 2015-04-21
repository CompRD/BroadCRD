///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#include "MainTools.h"
#include "Basevector.h"
#include "Qualvector.h"


static inline 
String Tag(String S = "MS") { return Date() + " (" + S + "): "; } 

int main(int argc, char *argv[])
{
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_LongLong(GENOME_SIZE);
  CommandArgument_UnsignedInt_OrDefault(COVERAGE, 50);
  CommandArgument_UnsignedInt_OrDefault(READ_LEN, 100);
  CommandArgument_String_OrDefault(HEAD_OUT, "");
  EndCommandArguments;
  
  {
    BaseVecVec reads;
    QualVecVec quals;
    
    const size_t n_reads = (GENOME_SIZE * COVERAGE) / READ_LEN;
    
    cout << 
    
    cout << Tag() << "Building a dummy set of " 
         << n_reads << " reads of length " << READ_LEN << " (expect 100 dots '.')." << endl;
    
    for (size_t i = 0; i != n_reads; i++) {
      if ((i + 1) * 100 / n_reads > i * 100 / n_reads) cout << "." << flush;
      
      BaseVec read(READ_LEN);
      QualVec qual(READ_LEN);
      
      for (size_t j = 0; j != READ_LEN; j++) {
        read.set(j, (i * (j + 7)) & 3u);  // silly formula
        qual[j] = (i * (j + 7)) % 40;  // silly formula
      }
      reads.push_back(read);
      quals.push_back(qual);
    }
    cout << endl;
    cout << Tag() << "Done building a dummy set of reads." << endl; 
    
    if (HEAD_OUT != "") {
      cout << Tag() << "Writing reads to disk." << endl; 
      reads.WriteAll(HEAD_OUT + ".fastb");
      cout << Tag() << "Writing quals to disk." << endl; 
      quals.WriteAll(HEAD_OUT + ".qualb");
    }
  }
  
  if (HEAD_OUT != "") {
    cout << Tag() << "Reading reads from disk." << endl; 
    BaseVecVec reads(HEAD_OUT + ".fastb");
    cout << Tag() << "Reading quals from disk." << endl; 
    QualVecVec quals(HEAD_OUT + ".qualb");
  }


  
  cout << Tag() << "Done." << endl; 

}

