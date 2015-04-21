/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// ------------------------------------------------------------------------
//
// KmerParcelsTest: A usage example for KmerParcels.
//
// 2009 Apt  Josh Burton
// 2010 Mar  Filipe Ribeiro
//
// ------------------------------------------------------------------------

#include "MainTools.h"
#include "Basevector.h"
#include "kmers/KmerParcels.h"



void LoadBases(const String & HEAD, const size_t n_reads, BaseVecVec * bases)
{
  if (HEAD != "") {
    bases->ReadAll(HEAD + ".fastb");
  }
  else {
    // Otherwise, build a sample BaseVecVec manually.
    // This will be a completely arbitrary dataset: 
    // 10000 reads, each with 40 bases, 
    // randomly generated each time KmerParcelsTest is run.
    for (size_t i = 0; i != n_reads; i++) {
      BaseVec b;
      b.SetFromStringWithNs("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");
      bases->push_back(b);
    }
    
    // Display this dataset.
    cout << "INPUT BASES: Here is our sample dataset." 
         << " It consists of " << n_reads << " reads, each with 40 randomly generated bases." 
         << endl << endl;
    for (size_t i = 0; i != n_reads; i++)
      (*bases)[i].Print(cout);
    cout << endl;
  }
}  






int main(int argc, char *argv[])
{
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_String_OrDefault_Doc(HEAD, "", "head of 1st .fastb file.  By default we randomly generate a BaseVecVec.");
  CommandArgument_UnsignedInt_OrDefault_Doc(K, 8, "Keep K=8 unless supplying your own FASTB.");
  CommandArgument_UnsignedInt_OrDefault(N_PARCELS, 8);
  CommandArgument_UnsignedInt_OrDefault(N_READS, 20);
  CommandArgument_Bool_OrDefault_Doc(WRITE_KMERS, False, "True: parcels on disk; False: parcels in memory." );
  EndCommandArguments;
  
  
  BaseVecVec bases; 
  
  LoadBases(HEAD, N_READS, &bases);
 

  // ---- create parcels
  cout << Date() << ": Splitting up this dataset into parcels..." << endl;


  KmerParcelsDiskStore parcels_disk(K, HEAD);
  KmerParcelsMemStore parcels_mem(K);

  KmerParcelsStore & parcels = (WRITE_KMERS) ? 
    static_cast<KmerParcelsStore &>(parcels_disk) : 
    static_cast<KmerParcelsStore &>(parcels_mem);
      
  
  const size_t n_threads = N_PARCELS;
  KmerParcelsBuilder builder(bases, parcels, n_threads);
  builder.Build();

  
  // ---- use parcels
  cout << Date() << ": Reading the kmer parcels back in..." << endl;

  for (size_t parcel_ID = 0; parcel_ID != N_PARCELS; parcel_ID++) {
    cout << Date() << ": Reading parcel # " << parcel_ID << endl;
    
    KmerParcelReader parcel_reader(parcels, parcel_ID);
    
    while (parcel_reader.GetNextKmerBatch()) {
      const KmerBatch & batch = parcel_reader.CurrentKmerBatch();
      
      batch.Print(cout);
      cout << endl;
    }
    cout << endl;
  }
  
  cout << Date() << ": Done!" << endl;
}
