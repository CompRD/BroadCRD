///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

//
//    Build unibases from kmer graph. 
//

#include "MainTools.h"
#include "Basevector.h"

#include "kmers/KmerSpectrumCore.h"

#include "kmers/naif_kmer/KmerFreqAffixesMap.h"
#include "kmers/naif_kmer/Unibases.h"


static inline 
String Tag(String S = "KTU") { return Date() + " (" + S + "): "; } 


bool bv_sz_lt(const BaseVec & bva, const BaseVec & bvb)
{
  if (bva.size() < bvb.size()) return true;
  if (bva.size() > bvb.size()) return false;
  return bva < bvb;
}


// ---- "Less Than" class for sorting BaseVec indexes

class IBVLT
{
  const BaseVecVec & _bvv;
public:
  IBVLT(const BaseVecVec & bvv) : _bvv(bvv) {}

  bool operator() (const size_t & ia, const size_t & ib) 
  {
    return bv_sz_lt(_bvv[ia], _bvv[ib]);
  }
};



// ---- Statistics class that stores, for each unibase, various metrics

class Stats 
{
public:
  vec<float> mean_kf;
  vec<unsigned> nb;
  vec<unsigned> nk;
  vec<bool> pal0;
  vec<bool> pal1;
  vec<bool> cycle;
  vec<unsigned> n_pre;
  vec<unsigned> n_suf;
  vec<unsigned> n_pre_suf;
  vec<unsigned> n_suf_pre;
  vec<unsigned> kf_pre_suf;
  vec<unsigned> kf_suf_pre;
  vec<size_t> iubs;

template<class KMER_REC_t>
  void analyze(const unsigned K, 
               const KmerMap<KMER_REC_t> & kmap,
               const BaseVecVec & bvv_ub)
  {
    typedef typename KMER_REC_t::kmer_type Kmer_t;

    const size_t nub = bvv_ub.size();
    for (size_t iub = 0; iub < nub; iub++) {
    
      const BaseVec & bv = bvv_ub[iub];
      const size_t nb = bv.size();
      const size_t nk = nb - K + 1;

      // ---- determine the number of prefixes of the first kmer 

      SubKmers<BaseVec, Kmer_t> kmer0_ub(K, bv, 0);
      const bool is_fw0 = kmer0_ub.is_canonical_fw();
      const bool is_pal0 = kmer0_ub.is_palindrome();
      const KMER_REC_t & krec0 = kmap(kmer0_ub.canonical());
      const unsigned n_pre = krec0.n_prefixes(is_fw0);

      // ---- determine the number of suffixes of the last kmer

      SubKmers<BaseVec, Kmer_t> kmer1_ub(K, bv, nb - K);
      const bool is_fw1 = kmer1_ub.is_canonical_fw();
      const bool is_pal1 = kmer1_ub.is_palindrome();
      const KMER_REC_t & krec1 = kmap(kmer1_ub.canonical());
      const unsigned n_suf = krec1.n_suffixes(is_fw1);

      // ---- determine the number of parallel kmer beginnings
      
      unsigned n_suf_pre = 0;
      unsigned kf_suf_pre = 0;
      KmerFWRC<Kmer_t> kmer0mFR = kmer0_ub;
      for (unsigned base = 0; base < 4; base++) {
        kmer0mFR.set(K - 1, base);
        const KMER_REC_t & kr = kmap(kmer0mFR.canonical());
        if (kr.is_valid_kmer()) {
          n_suf_pre++;
          kf_suf_pre += kr.freq();
        }
      }

      // ---- determine the number of parallel kmer ends
      
      unsigned n_pre_suf = 0;
      unsigned kf_pre_suf = 0;
      KmerFWRC<Kmer_t> kmer1mFR = kmer1_ub;
      for (unsigned base = 0; base < 4; base++) {
        kmer1mFR.set(0, base);
        const KMER_REC_t & kr = kmap(kmer1mFR.canonical());
        if (kr.is_valid_kmer()) {
          n_pre_suf++;
          kf_pre_suf += kr.freq();
        }
      }

      // ---- determine if unibase is a cycle or not

      bool is_cycle = false;
      for (unsigned i_suf = 0; i_suf < n_suf; i_suf++) {
        KmerFWRC<Kmer_t> kmer1sFR = kmer1_ub;
        kmer1sFR.push_right(krec1.suffix(i_suf, is_fw1));
        const KMER_REC_t & krec1s = kmap(kmer1sFR.canonical());
        if (krec1s.is_valid_kmer()) {
          if (krec1s == krec0)
            is_cycle = true;
        }
        else {
          cout << "oops! iub= " << iub << " i_suf=" << i_suf << "/" << n_suf << endl;
        }
        
      }


      this->nb.push_back(nb);
      this->nk.push_back(nk);
      this->pal0.push_back(is_pal0);
      this->pal1.push_back(is_pal1);
      this->cycle.push_back(is_cycle);
      this->n_pre.push_back(n_pre);
      this->n_suf.push_back(n_suf);
      this->n_pre_suf.push_back(n_pre_suf);
      this->n_suf_pre.push_back(n_suf_pre);
      this->kf_pre_suf.push_back(kf_pre_suf);
      this->kf_suf_pre.push_back(kf_suf_pre);
      this->iubs.push_back(iub);
    }
  }

  void sort_indexes_descending(const BaseVecVec & bvv)
  {
    IBVLT lt(bvv);
    sort(iubs.rbegin(), iubs.rend(), lt);
  }


  void to_text_file(const String & fn) 
  {
    ofstream outfs;
    outfs.open(fn.c_str());
    const String labels = ("#     i_ub"
                           "           nb"
                           "          cnb"
                           "           nk"
                           "          cnk"
                           "    mean_kf"
                           "   "
                           "  n_pre"
                           "  n_suf"
                           "   kf_suf_pre"
                           "   kf_pre_suf"
                           );
    outfs << labels << endl;
    
    const size_t nub = nb.size();
    size_t cnb = 0;
    size_t cnk = 0;
    for (size_t i = 0; i < nub; i++) {
      const size_t iub = iubs[i];
      cnb += nb[iub];
      cnk += nk[iub];
      outfs << setw(10) << i 
            << " " << setw(12) << nb[iub]
            << " " << setw(12) << cnb
            << " " << setw(12) << nk[iub]
            << " " << setw(12) << cnk
            << " " << setw(10) << fixed << setprecision(2) << mean_kf[iub]
            << "   " 
            << "  " << n_pre[iub] << " < " << n_suf_pre[iub]
            << "  " << n_pre_suf[iub] << " > " << n_suf[iub]
            << " " << setw(12) << kf_suf_pre[iub]
            << " " << setw(12) << kf_pre_suf[iub]
            << "  cycle=" << (cycle[iub] ? "yes" : "no ")
            << "  pal0="  << (pal0[iub]  ? "yes" : "no ")
            << "  pal1="  << (pal1[iub]  ? "yes" : "no ")
            << endl;
    }
    outfs << labels << endl;
    outfs.close();
  }

};






		     



template<class KMER_t>
void unibases_build(const unsigned K, 
                    const BaseVecVec & bvv, 
                    BaseVecVec * bvv_ub_p,
                    Stats * stats_ub_p,
                    KmerSpectrum *kspec_p,
                    KmerSpectrum *kspec_ub_p,
                    const Validator & validator,
                    const size_t verbosity,
                    const size_t n_threads)
{
  typedef KmerBVLoc<KmerFreqAffixes<KMER_t> > KmerRec_t;

  // ---- build the K-mer hash map

  cout << Tag() << "Building kmer map at K= " << K << "." << endl;

  KmerMap<KmerRec_t> kmap;
  const double hash_table_ratio = 1.5;
  kmer_freq_affixes_map_build_parallel(K, bvv, validator, hash_table_ratio, 
                                       & kmap,
                                       verbosity, n_threads, 0);
  
  // ---- verify kmap
  if (false) {
    cout << Tag() << "Verifying kmer map." << endl;
    cout << Tag() << "kmap.size()= " << kmap.size_hash() << endl;
    kmer_freq_affixes_map_verify(kmap);
  }
  
  // ---- build unibases
  
  cout << Tag() << "Building unibases at K= " << K << "." << endl;
  unibases_from_kmer_affix_map(K, & kmap, bvv_ub_p, &(stats_ub_p->mean_kf));


  // ---- kmer spectrum of reads from kmap

  if (kspec_p != 0) {
    kmer_spectrum_from_kmer_freq_map(kmap, kspec_p);
  }


  // ---- kmer spectrum of unibase kmers

  if (kspec_ub_p != 0) {
    const size_t nub = bvv_ub_p->size();
    for (size_t iub = 0; iub < nub; iub++) {
      float kf = stats_ub_p->mean_kf[iub];
      size_t kf0 = kf;
      size_t kf1 = kf + 1;
      size_t nk = (*bvv_ub_p)[iub].size() - K + 1;

      size_t n0 = 0.5 + float(nk) * (kf - float(kf0));
      size_t n1 = 0.5 + float(nk) * (float(kf1) - kf);
      if (kf1 >= kspec_ub_p->size()) kspec_ub_p->resize(kf1 + 1, 0);
      (*kspec_ub_p)[kf0] += n0;
      (*kspec_ub_p)[kf1] += n1;
    }
  }

  // ---- count kmers in unibases (just to make sure)

  const size_t nub = bvv_ub_p->size();
  size_t nk_tot = 0;
  for (size_t iub = 0; iub < nub; iub++) 
    nk_tot += (*bvv_ub_p)[iub].size() + 1 - K;
  
  cout << Tag() << setw(14) << nub << " unibases found." << endl;
  cout << Tag() << setw(14) << nk_tot << " total kmers in unibases." << endl;


  // ---- compute stats on each unibase

  cout << Tag() << setw(14) << "Computing unibase statistics." << endl;
  stats_ub_p->analyze(K, kmap, *bvv_ub_p);

  
  // ---- sort stuff
  
  cout << Tag() << setw(14) << "Sorting unibases by size." << endl;
  stats_ub_p->sort_indexes_descending(*bvv_ub_p);
  sort(bvv_ub_p->rbegin(), bvv_ub_p->rend(), bv_sz_lt);

}
















int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_UnsignedInt_Doc(K, "Kmer size.");
  CommandArgument_String_Doc(HEAD, "The base vectors filename.");
  CommandArgument_UnsignedInt_OrDefault_Doc(KF_MIN, 1, "Minimum kmer frequency.");
  CommandArgument_UnsignedInt_OrDefault_Doc(KF_MAX, 0, "Maximum kmer frequency.");
  CommandArgument_Bool_OrDefault_Doc(KF_MIN_FROM_SPECTRUM, False, "Use minimum from kmer spectrum for KF_MIN.");
  CommandArgument_Bool_OrDefault_Doc(DO_SPECTRA, False, "Whether to output kmer spectra.");
  

  CommandArgument_UnsignedInt_OrDefault_Doc(VERBOSITY, 1, "Verbosity level.");
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  EndCommandArguments;

  NUM_THREADS = configNumThreads(NUM_THREADS);
  
  if (K > 504) {
    cout << Tag() << "K= " << K << " too large. Not implemented." << endl;
    exit(1);
  }

  // ---- loading data
  
  const String fn_fastb = HEAD + ".fastb";
  cout << Tag() << "Loading base vectors from '" << fn_fastb << "'." << endl;
  BaseVecVec bvv(fn_fastb);
  KmerSpectrum kspec(K);
  KmerSpectrum kspec_ub(K);

  if (KF_MIN_FROM_SPECTRUM) {
    cout << Tag() << "Computing kmer spectrum of reads." << endl;
    kmer_spectrum_compute(bvv, &kspec, VERBOSITY, NUM_THREADS);
    const size_t nbv = bvv.size();
    size_t read_len = 0;
    for (size_t ibv = 0; ibv < nbv; ibv++)
      read_len += bvv[ibv].size();
    read_len /= nbv;
    const unsigned ploidy = 2;
    const unsigned kf_low = 3;
    genome_analysis_report(kspec, read_len, ploidy, kf_low, VERBOSITY);
    KF_MIN = kspec.kf_min1();
    kspec.to_text_file(HEAD + ".ufk");
  }

  // ---- building unibases

  cout << Tag() << "Computing unibases at K= " << K << "." << endl;
  BaseVecVec bvv_ub;
  Stats stats_ub;
  
  Validator val(KF_MIN, KF_MAX);
  cout << Tag() << "KF_MIN= " << KF_MIN << endl;
  cout << Tag() << "KF_MAX= " << KF_MAX << endl;


  KmerSpectrum * kspec_p    = (DO_SPECTRA) ? & kspec : 0;
  KmerSpectrum * kspec_ub_p = (DO_SPECTRA) ? & kspec_ub : 0;
  
  if      (K <=  29) unibases_build<Kmer29> (K, bvv, & bvv_ub, & stats_ub, kspec_p, kspec_ub_p, val, VERBOSITY, NUM_THREADS);
  else if (K <=  60) unibases_build<Kmer60> (K, bvv, & bvv_ub, & stats_ub, kspec_p, kspec_ub_p, val, VERBOSITY, NUM_THREADS);
  else if (K <= 124) unibases_build<Kmer124>(K, bvv, & bvv_ub, & stats_ub, kspec_p, kspec_ub_p, val, VERBOSITY, NUM_THREADS);
  else if (K <= 248) unibases_build<Kmer248>(K, bvv, & bvv_ub, & stats_ub, kspec_p, kspec_ub_p, val, VERBOSITY, NUM_THREADS);
  else if (K <= 504) unibases_build<Kmer504>(K, bvv, & bvv_ub, & stats_ub, kspec_p, kspec_ub_p, val, VERBOSITY, NUM_THREADS);


  // ---- save kmer spectra

  if (DO_SPECTRA) {
    kspec.to_text_file(HEAD + ".ufk");
    kspec_ub.to_text_file(HEAD + ".ufk.ub");
  }

  // ---- save unibases

  const String fn_ub_fastb = HEAD + "." + ToString(K) + "mer.ub.fastb";
  
  cout << Tag() << "Saving unibases to '" << fn_ub_fastb << "'." << endl;
  bvv_ub.WriteAll(fn_ub_fastb);

  
  const String fn_ub_stats = HEAD + "." + ToString(K) + "mer.ub.stats";

  cout << Tag() << "Saving unibase stats to '" << fn_ub_stats << "'." << endl;
  stats_ub.to_text_file(fn_ub_stats);

}


