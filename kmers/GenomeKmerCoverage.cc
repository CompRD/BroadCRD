///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#include "MainTools.h"

#include "Basevector.h"

static inline 
String Tag(String S = "GKC") { return Date() + " (" + S + "): "; } 

#include "kmers/naif_kmer/NaifKmerizer.h" 
#include "kmers/naif_kmer/KernelKmerStorer.h" 
#include "kmers/naif_kmer/KmerMap.h" 






template<class MAP_t>
void genome_kmer_coverage_analyze(const size_t K,
                                  const BaseVecVec & bases,
                                  const MAP_t & kmap,
                                  const String & fn)
{
  typedef map<size_t, vec<size_t> > UBSizesMap;

  UBSizesMap ub_sizes_of_cov;

  bool to_file = false;
  ofstream os;
  if (to_file) os.open(fn.c_str());

  KmerKmerFreq<Kmer29> kmerkf(0);

  const size_t nbv = bases.size();
  for (size_t ibv = 0; ibv != nbv; ibv++) {
    size_t c_prev = 0;
    size_t ub_sz = 0;
    SubKmers<BaseVec, Kmer29> kmer_cur(K, bases[ibv]);
    while (kmer_cur.not_done()) {
      Kmer29 & kmer = kmer_cur.canonical();

      kmerkf = kmap(kmer);
      if (kmerkf.is_valid_kmer()) {

	if (to_file) os << kmerkf.freq() << endl;
	
	if (kmerkf.freq() == c_prev) 
	  ub_sz++;
	else {
	  if (ub_sz > 0) 
	    ub_sizes_of_cov[c_prev].push_back(ub_sz);
	  ub_sz = 1;
	  c_prev = kmerkf.freq();
	}
      }
      kmer_cur.next();
    }
    if (to_file) os << endl << endl;
  }

  /*

  for (size_t i = 0; i != nr; i++) {
    const BaseVec & bv = bases[i];
    const size_t nb = bv.size();

    if (nb >= K) {

      Kmer kmerFW(K);
      Kmer kmerRC(K);
      
      // ---- load the 1st K - 1 bases
      for (unsigned ib = 0; ib != K - 1; ib++) {
	const uint64_t base = bv[ib];
	kmerFW.push_right(base);
	kmerRC.push_left(3ul ^ base);  // the complement of base
      }
      
      
      size_t c_prev = 0;
      size_t ub_sz = 0;

      // ---- build kmers by adding one base at a time
      for (unsigned ib = K - 1; ib != nb; ib++) {
	const uint64_t base = bv[ib];
	kmerFW.push_right(base);
	kmerRC.push_left(3ul ^ base);  // the complement of base
	
	const Kmer & kmer = (kmerFW < kmerRC) ? kmerFW : kmerRC; 
	
        kmap.binary_search(kmer, &kkf);
        if (to_file) os << kkf.freq() << endl;

        if (kkf.freq() == c_prev) 
          ub_sz++;
        else {
          if (ub_sz > 0) 
            ub_sizes_of_cov[c_prev].push_back(ub_sz);
          ub_sz = 1;
          c_prev = kkf.freq();
        }

      }
      if (to_file) os << endl << endl;
    }
  }
  */
  for (UBSizesMap::iterator it = ub_sizes_of_cov.begin();
       it != ub_sizes_of_cov.end(); it++) {
    vec<size_t> & sizes = it->second;
    sort(sizes.begin(), sizes.end());
    size_t n = sizes.size();
    size_t nb = 0;
    for (size_t i = 0; i != n; i++)
      nb += sizes[i];
    
    size_t cum = sizes[0];
    size_t iN50 = 0;
    while (cum < nb / 2) {
      iN50++;
      cum += sizes[iN50];
    }
    size_t N50 = sizes[iN50];
      
    size_t cov = it->first;
    cout << "cov= " << setw(7)  << cov << "   "
         << "n= "   << setw(7)  << n   << "   "
         << "nb= "  << setw(10) << nb  << "   "
         << "n50= " << N50
         << endl;
  }



  if (to_file) os.close();

}








int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_UnsignedInt_Doc(K, "Kmer size.");
  CommandArgument_String_Doc(HEAD, "looks for reference genome in <HEAD>.fastb");
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, "Number of threads to use in parallel");

  EndCommandArguments;

  NUM_THREADS = configNumThreads(NUM_THREADS);
 

  
  cout << Tag() << "Loading bases" << endl;
  BaseVecVec bases(HEAD + ".fastb");

  const bool verbose = false;

  cout << Tag() << "Building Kmer vector" << endl;
  vec<KmerKmerFreq<Kmer29> > kvec;
  {
    KernelKmerStorer<KmerKmerFreq<Kmer29> > storer(bases, K, &kvec);
    naif_kmerize(&storer, NUM_THREADS, verbose);
  }

  cout << Tag() << "Building Kmer map" << endl;

  sort(kvec.begin(), kvec.end(), &(kmer_freq_gt<KmerKmerFreq<Kmer29> >));
  KmerMap<KmerKmerFreq<Kmer29> > kmap(kvec);
  //KmerMapBin<KmerKmerFreq> kmap(kvec);

  //kmap.output_hash_stats();

  
  cout << Tag() << "Processing genome" << endl;

  genome_kmer_coverage_analyze(K, bases, kmap, HEAD + ".K=" + ToString(K) + ".kc");

  cout << Tag() << "Done" << endl;

}
