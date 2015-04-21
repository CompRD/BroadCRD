///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#include <stdlib.h>   // rand() in [0, RAND_MAX]

#include "MainTools.h"


#include "Basevector.h"
#include "PairsManager.h"
#include "Qualvector.h"
#include "feudal/IncrementalWriter.h"

static inline 
String Tag(String S = "RFBKF") { return Date() + " (" + S + "): "; } 


#include "kmers/naif_kmer/NaifKmerizer.h" 
#include "kmers/naif_kmer/KernelReadKmerFreqFinder.h" 






int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String_Doc(HEAD, "looks for <HEAD>.fastb");
  CommandArgument_String_OrDefault_Doc(HEAD_OUT, "", "writes kept read set to <HEAD_OUT>.fastb");
  CommandArgument_String_OrDefault_Doc(HEAD_OUT2, "", "writes unkept reads to <HEAD_OUT2>.fastb");
  CommandArgument_UnsignedInt_OrDefault_Doc(K, 27, "Kmer size.");
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, "Number of threads to use in parallel");
  CommandArgument_UnsignedInt_Doc(KF_MAX, "Kmer frequency maximum above which to downsample reads.");
  CommandArgument_LongLong_OrDefault(MAX_MEMORY_GB, 0);
  EndCommandArguments;

  SetMaxMemory(MAX_MEMORY_GB << 30);

  if (true) {
    const size_t n_cpus = processorsOnline();
    if (NUM_THREADS == 0 || NUM_THREADS > n_cpus)
      NUM_THREADS = n_cpus;
  }

  bool verbose = false;

  String head_out = (HEAD_OUT == "") ? (HEAD + ".kflt_" + ToString(KF_MAX)) : HEAD_OUT;

  cout << Tag() << "Loading reads, quals, & pairs" << endl;
  PairsManager pairs(HEAD + ".pairs");
  const BaseVecVec bases_in(HEAD + ".fastb");
  VirtualMasterVec<QualVec> quals_in((HEAD + ".qualb").c_str());

  const size_t n_bv = bases_in.size();

  // ---- Compute max kmer freq per read
  vec<KFStats> kf_stats(n_bv);
  ReadKmerFreqFinder<Kmer29> finder(K, bases_in, &kf_stats, NUM_THREADS);
  naif_kmerize(&finder, NUM_THREADS, verbose);

  // ---- Generate output
  
  IncrementalWriter<BaseVec> bases_out((head_out + ".fastb").c_str());
  IncrementalWriter<QualVec> quals_out((head_out + ".qualb").c_str());
  PairsManager pairs_out;

  cout << Tag() << setw(12) << n_bv << " reads in" << endl;
  size_t n_bv_out = 0;
  vec<bool> pairs_kept(pairs.nPairs(), false);
  for (size_t p = 0; p < pairs.nPairs(); p++) {
    size_t i1 = pairs.ID1(p), i2 = pairs.ID2(p);
    bool add = true;
    const size_t kf_max_pair = Min(kf_stats[i1].max, kf_stats[i2].max);
    if (kf_max_pair > KF_MAX) {
      add = ((rand() % kf_max_pair) < KF_MAX) ? true : false;
    }

    if (add) {
      bases_out.add(bases_in[i1]);
      bases_out.add(bases_in[i2]);
      quals_out.add(quals_in[i1]);
      quals_out.add(quals_in[i2]);
      pairs_out.addPair(n_bv_out, n_bv_out + 1,
			pairs.sep(p), pairs.sd(p),
			pairs.libraryName(p), True);
      n_bv_out += 2;
      pairs_kept[p] = true;
    }
  }
  bases_out.close();
  quals_out.close();
  pairs_out.Write(head_out + ".pairs");

  cout << Tag() << setw(12) << n_bv_out << " reads written to " << head_out << endl;

  int n_bv_out2 = 0;
  if (HEAD_OUT2 != "") {
    IncrementalWriter<BaseVec> bases_out2((HEAD_OUT2 + ".fastb").c_str());
    IncrementalWriter<QualVec> quals_out2((HEAD_OUT2 + ".qualb").c_str());
    PairsManager pairs_out2;
    for (size_t p = 0; p < pairs.nPairs(); p++) {
      if (!pairs_kept[p]) {
	size_t i1 = pairs.ID1(p), i2 = pairs.ID2(p);
	bases_out2.add(bases_in[i1]);
	bases_out2.add(bases_in[i2]);
	quals_out2.add(quals_in[i1]);
	quals_out2.add(quals_in[i2]);
	pairs_out2.addPair(n_bv_out2, n_bv_out2 + 1,
			  pairs.sep(p), pairs.sd(p),
			  pairs.libraryName(p), True);
	n_bv_out2 += 2;
      }
    }
    pairs_out2.Write(HEAD_OUT2 + ".pairs");
    cout << Tag() << setw(12) << n_bv_out2 << " reads written to " << HEAD_OUT2 << endl;
  }
}
