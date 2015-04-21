///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/// Program: CompareKmerParcels
///
/// Takes two sets of kmer parcels and compares them
///


#include "Basevector.h"
#include "MainTools.h"
#include "kmers/KmerParcels.h" 
#include "kmers/KmerParcelsTools.h" 

#include <map>


static inline 
String Tag(String S = "KPC") { return Date() + " (" + S + "): "; } 



int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_UnsignedInt_Doc(K, "Kmer size");
  CommandArgument_String_Doc(HEAD1, "looks for directory <HEAD1>.<K>merParcels");
  CommandArgument_String_Doc(HEAD2, "looks for directory <HEAD2>.<K>merParcels");
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 4, "number of parallel threads to use.");
  EndCommandArguments;


  cout << Tag() << "BEGIN" << endl;

  KmerParcelsDiskStore parcels1(K, HEAD1);
  KmerParcelsDiskStore parcels2(K, HEAD2);
  
  CompareKmerParcels(parcels1, parcels2, HEAD1, HEAD2, NUM_THREADS);
                       
  cout << Tag() << "END" << endl;
}
