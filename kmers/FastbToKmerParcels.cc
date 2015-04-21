/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/// Program: FastbToKmerParcels
/// 
///  parcelizes a .fastb file in kmer parcels


#include "Basevector.h"
#include "MainTools.h"
#include "kmers/KmerParcels.h" 

static inline 
String Tag(String S = "Fb2KP") { return Date() + " (" + S + "): "; } 

int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_Int_Doc
    (K, "Size of Kmer for KmerParcel building");
  CommandArgument_String_Doc
    (HEAD, "looks for HEAD.fastb, creates HEAD.<K>merParcels/<parcels>");
  CommandArgument_Int_OrDefault_Doc
    (NUM_THREADS, 8, "Number of threads to use in parallel" );
  CommandArgument_Int_OrDefault_Doc
    (NUM_PARCELS, 0, "Number of parcels" );
  CommandArgument_Int_OrDefault_Doc
    (MIN_KMER_FREQ, 1, "The minimum kmer frequency (or batch size) to include in parcels" );
  CommandArgument_Bool_OrDefault_Doc
    (WRITE_PARCELS, True, "Write parcels to disk.  Useful if you DON'T want the parcels written to disk, just the statistics." );
  CommandArgument_Bool_OrDefault_Doc
    (WRITE_GC_STATS, True, "Write various kmer frequency counts based on kmer GC content." );
  CommandArgument_Bool_OrDefault_Doc
    (CANONICALIZE, True, "To canonicalize or not to canonicalize, that is the question.  Useful if you want to look at coverage assymmetry." );
  EndCommandArguments;
  
  cout << Tag() << "START: FastbToKmerParcels" << endl;

  cout << Tag() << "Loading bases" << endl;
  vecbasevector bases(HEAD + ".fastb");
  
  if (bases.empty())
    FatalErr("No bases to parcel - the fastb file is empty!");

  cout << Tag() << "Building kmer parcels" << endl;
  

  KmerParcelsDiskStore parcels(K, HEAD);

  if (WRITE_PARCELS) { 
    // the destructor of parcels will not erase them from disk 
    parcels.SetKeepOnDisk(true);
  }
  else { // we only want statistics
    // the destructor of parcels will erase them from disk 
    parcels.SetKeepOnDisk(false);
    // the kmer validator makes sure that no kmers are written to disk
    parcels.KmerValidator().SetNone();
  }

  


  KmerParcelsBuilder parcel_builder(bases, parcels, NUM_THREADS);
  
  parcel_builder.SetVerbose(true);
  parcel_builder.SetCanonicalize(CANONICALIZE);
  parcel_builder.SetWriteGCStats(WRITE_GC_STATS);

  parcel_builder.Build(NUM_PARCELS); // NUM_PARCELS = 0 => automatically compute num_parcels
  size_t num_parcels = parcel_builder.GetNumParcels();
  

  cout << Tag() << "GREPABLE: num_parcels = " << num_parcels << endl;
  
  {
    GenomeSizeEstimator gse(parcel_builder.GetKmerFrequencyCounters());
    
    const size_t G  = gse.genome_size_total();
    const size_t G1 = gse.genome_size_unique();
    const size_t GR = gse.genome_size_repetitive();
    
    cout << Tag() << "GREPABLE: Genome size estimate        = " 
         << setw(14) << ToStringAddCommas(G) << " bases" << endl;
    cout << Tag() << "GREPABLE: Genome size estimate CN = 1 = " 
         << setw(14) << ToStringAddCommas(G1) << " bases ( " 
         << fixed << setw(5) << setprecision(1) << 100 * float(G1)/float(G) << " % )" << endl;
    cout << Tag() << "GREPABLE: Genome size estimate CN > 1 = " 
         << setw(14) << ToStringAddCommas(GR) << " bases ( "
         << fixed << setw(5) << setprecision(1) << 100 * float(GR)/float(G) << " % )" << endl;
    
    const size_t nr = bases.size();
    size_t nb = 0;
    for (size_t i = 0; i != nr; i++)
      nb += bases[i].size();
    const float coverage = float(nb) / float(G);
    cout << Tag() << "GREPABLE: Coverage estimate = " << setprecision(1) << coverage << " x" << endl;
  }


  cout << Tag() << "END: FastbToKmerParcels" << endl;
}
