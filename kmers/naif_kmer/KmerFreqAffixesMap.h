///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#ifndef KMERS__NAIF_KMER__KMER_FREQ_AFFIXES_MAP_H
#define KMERS__NAIF_KMER__KMER_FREQ_AFFIXES_MAP_H


#include "kmers/KmerSpectra.h"

#include "kmers/naif_kmer/NaifKmerizer.h" 
#include "kmers/naif_kmer/KernelKmerStorer.h" 
#include "kmers/naif_kmer/KmerMap.h" 
#include "system/WorklistN.h"


static inline 
String TagKFAM(String S = "KFAM") { return Date() + " (" + S + "): "; } 

class FreqAffixes
{
  uint64_t _freq : 56;
  uint64_t _pre  :  4;  // 1 bit per possible base
  uint64_t _suf  :  4;  // 1 bit per possible base

  unsigned _n_bases(const unsigned bit4) const 
  { 
    return ((bit4 & 1u) + ((bit4 >> 1) & 1u) + ((bit4 >> 2) & 1u) + ((bit4 >> 3) & 1u)); 
  }

  unsigned _base(const unsigned bit4, const unsigned i) const
  {
    unsigned n = 0;
    if (bit4 & 1u) n++;   if (n > i) return 0u;
    if (bit4 & 2u) n++;   if (n > i) return 1u;
    if (bit4 & 4u) n++;   if (n > i) return 2u;
    if (bit4 & 8u) n++;   if (n > i) return 3u;
    return 4u;
  }

  String _str_affixes(const unsigned aff, const bool fw) const 
  {
    String s = "[" + (((aff & 1u) ? hieroglyph(fw ? 0 : 3) : " ") +
		      ((aff & 2u) ? hieroglyph(fw ? 1 : 2) : " ") +
		      ((aff & 4u) ? hieroglyph(fw ? 2 : 1) : " ") +
		      ((aff & 8u) ? hieroglyph(fw ? 3 : 0) : " ")) + "]";
    return s;
  }



public:
  FreqAffixes() : _freq(0), _pre(0), _suf(0) {}

  void set_freq  (const size_t freq)    { _freq = freq; }
  void set_prefix(const unsigned base)  { _pre |= (1u << (base & 3)); }
  void set_suffix(const unsigned base)  { _suf |= (1u << (base & 3)); }


  size_t freq() const { return _freq; }

  bool has_prefix(const unsigned base) const { return _pre & (1u << (base & 3)); }
  bool has_suffix(const unsigned base) const { return _suf & (1u << (base & 3)); }

  size_t n_prefixes(const bool fw = true) const 
  { return (fw ? _n_bases(_pre) : _n_bases(_suf)); } 

  size_t n_suffixes(const bool fw = true) const
  { return (fw ? _n_bases(_suf) : _n_bases(_pre)); }

  size_t n_affixes() const 
  { return _n_bases(_pre) + _n_bases(_suf); }


  unsigned prefix(const unsigned i, const bool fw = true) const
  {
    return (fw ? _base(_pre, i) : 3u ^ _base(_suf, i));
  }

  unsigned suffix(const unsigned i, const bool fw = true) const
  {
    return (fw ? _base(_suf, i) : 3u ^ _base(_pre, i));
  }

  String str_prefixes(const bool fw = true) const 
  {
    return _str_affixes(fw ? _pre : _suf, fw);
  }

  String str_suffixes(const bool fw = true) const 
  {
    return _str_affixes(fw ? _suf : _pre, fw);
  }
};





template<class KMER_t>
class KmerFreqAffixes : public KMER_t, public FreqAffixes
{
public:
  KmerFreqAffixes(const KMER_t & kmer) : KMER_t(kmer), FreqAffixes() {}
  explicit KmerFreqAffixes(const unsigned K = 0) : KMER_t(K), FreqAffixes() {}
 
  friend
  bool operator < (const KmerFreqAffixes & a, const KmerFreqAffixes & b)
  { 
    return (static_cast<const KMER_t &>(a) < static_cast<const KMER_t &>(b));
  }
 
};





TRIVIALLY_SERIALIZABLE(KmerFreqAffixes<Kmer29>);
TRIVIALLY_SERIALIZABLE(KmerFreqAffixes<Kmer60>);
TRIVIALLY_SERIALIZABLE(KmerFreqAffixes<Kmer124>);
TRIVIALLY_SERIALIZABLE(KmerFreqAffixes<Kmer248>);
TRIVIALLY_SERIALIZABLE(KmerFreqAffixes<Kmer504>);













template<class KMER_REC_t>
class KmerAffixesMapNavigator
{
  typedef typename KMER_REC_t::kmer_type Kmer_t;

  const KmerMap<KMER_REC_t> & _kmap;
  Kmer_t _kmer_fw;
  Kmer_t _kmer_rc;
  bool _fw;
  KMER_REC_t _kmer_rec;

public:
  KmerAffixesMapNavigator(const KmerMap<KMER_REC_t> & kmap,
			  const Kmer_t & kmer) : 
    _kmap(kmap),
    _kmer_fw(kmer),
    _kmer_rc(reverse_complement(kmer)),
    _fw(_kmer_fw < _kmer_rc),
    _kmer_rec(kmap(_fw ? _kmer_fw : _kmer_rc))
  {
    ForceAssert(_kmer_rec.is_valid_kmer());
  }

  KMER_REC_t & rec() const  { return _kmer_rec; } 
  
  Kmer_t & fw() const { return _kmer_fw; }
  Kmer_t & rc() const { return _kmer_rc; }


  unsigned n_prefixes() const { return _kmer_rec.n_prefixes(_fw); }
  unsigned n_suffixes() const { return _kmer_rec.n_suffixes(_fw); }

  unsigned prefix(unsigned i_pre) const { return _kmer_rec.prefix(i_pre, _fw); }
  unsigned suffix(unsigned i_suf) const { return _kmer_rec.suffix(i_suf, _fw); }

  void next_by_prefix(const unsigned i_pre)
  {
    ForceAssertLt(i_pre, n_prefixes());

    unsigned base_fw = prefix(i_pre);
    unsigned base_rc = 3u ^ base_fw;
    _kmer_fw.push_left (base_fw);
    _kmer_rc.push_right(base_rc);
    _fw = (_kmer_fw < _kmer_rc);
    _kmer_rec = _kmap(_fw ? _kmer_fw : _kmer_rc);

    ForceAssert(_kmer_rec.is_valid_kmer());
  }

  void next_by_suffix(const unsigned i_suf)
  {
    ForceAssertLt(i_suf, n_suffixes());

    unsigned base_fw = suffix(i_suf);
    unsigned base_rc = 3u ^ base_fw;
    _kmer_fw.push_right(base_fw);
    _kmer_rc.push_left (base_rc);
    _fw = (_kmer_fw < _kmer_rc);
    _kmer_rec = _kmap(_fw ? _kmer_fw : _kmer_rc);

    ForceAssert(_kmer_rec.is_valid_kmer());
  }
};








/*
template<class KMER_t>
void navigate_suffixes_up_to(const KMER_t kmer_final_fw,
                             BaseVec * bv_p,
                             const size_t nb_max,
                             const size_t n_branches_max) 
{
  if (kmer_final_fw != _kmer_fw &&
      n_branches_max > 0 &&
      bv_p->size() < nb_max) {
    
    while (n_suffixes() == 1 && 
           kmer_final_fw != _kmer_fw &&
           bv_p->size() < nb_max) {
      
      
      bv_p.push_back(suffix(0));
      next_by_suffix(0);
      
      
      
    }
    
    const unsigned n_suf = n_suffixes();
    for (unsigned i_suf = 0; i_suf < n_suf; i_suf++) {
      BaseVec bv;
      
      
      while (bv.size() < nb_max) {
        if (n_suffixes() == 0);
      }
    }
    
  }
}

*/











template<class KMER_REC_t>
class AffixesAddProc
{
  typedef typename KMER_REC_t::kmer_type Kmer_t;
  
  KmerMap<KMER_REC_t> & _kmap;
  const size_t _n_threads;
  const unsigned _verbosity;

public:
  AffixesAddProc(KmerMap<KMER_REC_t> * kmap_p, 
                 const size_t n_threads,
                 const unsigned verbosity) :
    _kmap(*kmap_p),
    _n_threads(n_threads),
    _verbosity(verbosity)
  {}

  AffixesAddProc(const AffixesAddProc & that) :
    _kmap(that._kmap),
    _n_threads(that._n_threads),
    _verbosity(that._verbosity)
  {}

  void operator() (const size_t i_thread)
  {
    const size_t nh = _kmap.size_hash();
    const size_t ih0 =  i_thread      * nh / _n_threads;
    const size_t ih1 = (i_thread + 1) * nh / _n_threads;

    //cout << "nt= " << _n_threads << " it= " << i_thread << " ih0= " << ih0 << " ih1= " << ih1 << endl;

    for (size_t ih = ih0; ih < ih1; ih++) {
      if (i_thread == 0 && _verbosity > 0)
        dots_pct(ih, ih1);

      KMER_REC_t & krec0 = _kmap[ih];
      if (krec0.is_valid_kmer()) { 
        
        // ---- search for prefixes
        {
          KmerFWRC<Kmer_t> kmerFR(krec0);
          kmerFR.push_left(0);
          
          for (unsigned base = 0; base < 4; base++) {
            kmerFR.set(0, base);
            if (_kmap(kmerFR.canonical()).is_valid_kmer())
              krec0.set_prefix(base);
          }
        }

        // ---- search for suffixes
        {
          const unsigned K = krec0.K();
          KmerFWRC<Kmer_t> kmerFR(krec0);
          kmerFR.push_right(0);
          
          for (unsigned base = 0; base < 4; base++) {
            kmerFR.set(K - 1, base);
            if (_kmap(kmerFR.canonical()).is_valid_kmer())
              krec0.set_suffix(base);
          }
        }
        
      }
    }      
    
  }

};






template<class KMER_REC_t>
void kmer_freq_affixes_map_build_parallel(const size_t K,
                                          const BaseVecVec & bvv,
                                          const Validator & validator_kf,
                                          const double hash_table_ratio,
                                          KmerMap<KMER_REC_t> * kmap_p,
                                          const size_t verbosity,
                                          const size_t n_threads,
                                          const size_t mean_mem_ceil = 0)
{
  const bool do_affixes = false;

  // ---- build kmer vector with frequencies and affixes
  
  vec<KMER_REC_t> kvec;
  
  if (do_affixes) {
    KernelKmerAffixesStorer<KMER_REC_t> storer(bvv, K, &kvec, &validator_kf);
    naif_kmerize(&storer, n_threads, verbosity, mean_mem_ceil);
  }
  else {
    KernelKmerStorer<KMER_REC_t> storer(bvv, K, &kvec, &validator_kf);
    naif_kmerize(&storer, n_threads, verbosity, mean_mem_ceil);
  }


  if (verbosity > 0)
    cout << TagKFAM() << setw(14) << kvec.size() << " records found." << endl;
  
  
  // ---- sort kvec by highest kmer frequency 
  //      since we are building a chain hash we want to add first the 
  //      high frequency kmers so that, when we recall them (which will happen
  //      often) they'll come up first.
  
  if (verbosity > 0) cout << TagKFAM() << "Sorting records." << endl; 
  sort(kvec.begin(), kvec.end(), &(kmer_freq_gt<KMER_REC_t>));

  // ---- convert from vec<kmer> to hash table
  
  if (verbosity > 0) cout << TagKFAM() << "Building " << K << "-mer hash table." << endl;
  kmap_p->from_kmer_vec(kvec, hash_table_ratio, verbosity);
  

  if (!do_affixes) {
    if (verbosity > 0) cout << TagKFAM() << "Finding affixes in parallel." << endl;
    
    AffixesAddProc<KMER_REC_t> adder(kmap_p, n_threads, verbosity);
    if (n_threads <= 1) adder(0);
    else
      parallelFor(0ul,n_threads,adder,n_threads);
  }


  if (verbosity > 0) cout << TagKFAM() << "Done building " << K << "-mer hash table." << endl;
  
  if (verbosity > 0) 
    kmer_affixes_map_freq_table_print(*kmap_p);

  //kmer_affixes_map_verify(*kmap_p);

}








template <class KMER_REC_t>
void kmer_freq_affixes_map_build_parallel(const size_t K,
                                          const BaseVecVec & bvv,
                                          const double hash_table_ratio,
                                          KmerMap<KMER_REC_t> * kmap_p,
                                          const size_t verbosity,
                                          const size_t n_threads,
                                          const size_t mean_mem_ceil = 0)
{
  Validator validator_kf;
  kmer_freq_affixes_map_build_parallel(K, bvv, validator_kf, hash_table_ratio, 
                                       kmap_p,
                                       verbosity, n_threads, mean_mem_ceil);
}








template <class KMER_REC_t>
void kmer_spectrum_from_kmer_freq_map(const KmerMap<KMER_REC_t> & kmap,
                                      KmerSpectrum * kspec_p)
{
  const size_t nh = kmap.size_hash();
  for (size_t ih = 0; ih != nh; ih++) {
    const KMER_REC_t & krec = kmap[ih];
    if (krec.is_valid_kmer()) {
      const size_t freq = krec.freq();
      if (kspec_p->size() <= freq)
        kspec_p->resize(freq + 1, 0);
      (*kspec_p)[freq]++;
    }
  }
}

  



// ---- print a table of the kmer frequency regarding number of affixes

template<class KMER_REC_t>
void kmer_affixes_map_freq_table_print(const KmerMap<KMER_REC_t> & kmap)
{
  const size_t nh = kmap.size_hash();
  vec<vec<size_t> > n_kmers(5, vec<size_t>(5, 0));
  vec<vec<size_t> > kf(5, vec<size_t>(5, 0));
    
  size_t n_kmers_total = 0;

  for (size_t ih = 0; ih < nh; ih++) {
    const KMER_REC_t & krec = kmap[ih];
    if (krec.is_valid_kmer()) { 
      unsigned n_pre = krec.n_prefixes();
      unsigned n_suf = krec.n_suffixes();
      n_kmers[n_pre][n_suf]++;
      n_kmers_total++;
      kf[n_pre][n_suf] += krec.freq();
    }
  }
  // ---- output table of flows
  cout << TagKFAM() << endl;
    
  for (size_t n_pre = 0; n_pre <= 4; n_pre++) {
    cout << TagKFAM() << "n_kmers(" << n_pre << "-" << n_pre << ")=  " 
	 << setw(10) << n_kmers[n_pre][n_pre]
	 << " (mean_freq= " << setw(6) << kf[n_pre][n_pre] / (n_kmers[n_pre][n_pre] + 1) << ")"
	 << endl;
    for (size_t n_suf = n_pre + 1; n_suf <= 4; n_suf++) {
      size_t nk = (n_kmers[n_pre][n_suf] + n_kmers[n_suf][n_pre]);
      size_t kf2 = (kf[n_pre][n_suf] + kf[n_suf][n_pre]);
      cout << TagKFAM() << "n_kmers(" << n_pre << "-" << n_suf << ")=  " 
	   << setw(10) << nk
	   << " (mean_freq= " << setw(6) << kf2 / (nk + 1) << ")"
	   << endl;
    }
    cout << TagKFAM() << endl;
  }
  cout << TagKFAM() << "n_kmers(tot)=  " << setw(10) << n_kmers_total << endl;
  cout << TagKFAM() << endl;
}





template<class KMER_REC_t>
void kmer_affixes_vec_freq_table_print(const vec<KMER_REC_t> & kvec)
{
  const size_t nh = kvec.size();
  vec<vec<size_t> > n_kmers(5, vec<size_t>(5, 0));
  vec<vec<size_t> > kf(5, vec<size_t>(5, 0));
    
  size_t n_kmers_total = 0;

  for (size_t ih = 0; ih < nh; ih++) {
    const KMER_REC_t & krec = kvec[ih];
    unsigned n_pre = krec.n_prefixes();
    unsigned n_suf = krec.n_suffixes();
    n_kmers[n_pre][n_suf]++;
    n_kmers_total++;
    kf[n_pre][n_suf] += krec.freq();
 }
  // ---- output table of flows
  cout << TagKFAM() << endl;
    
  for (size_t n_pre = 0; n_pre <= 4; n_pre++) {
    cout << TagKFAM() << "n_kmers(" << n_pre << "-" << n_pre << ")=  " 
	 << setw(10) << n_kmers[n_pre][n_pre]
	 << " (mean_freq= " << setw(6) << kf[n_pre][n_pre] / (n_kmers[n_pre][n_pre] + 1) << ")"
	 << endl;
    for (size_t n_suf = n_pre + 1; n_suf <= 4; n_suf++) {
      size_t nk = (n_kmers[n_pre][n_suf] + n_kmers[n_suf][n_pre]);
      size_t kf2 = (kf[n_pre][n_suf] + kf[n_suf][n_pre]);
      cout << TagKFAM() << "n_kmers(" << n_pre << "-" << n_suf << ")=  " 
	   << setw(10) << nk
	   << " (mean_freq= " << setw(6) << kf2 / (nk + 1) << ")"
	   << endl;
    }
    cout << TagKFAM() << endl;
  }
  cout << TagKFAM() << "n_kmers(tot)=  " << setw(10) << n_kmers_total << endl;
  cout << TagKFAM() << endl;
}








// ---- Makes sure that the affixes make sense



template<class KMER_REC_t>
void kmer_freq_affixes_map_verify(const KmerMap<KMER_REC_t> & kmap)
{
  typedef typename KMER_REC_t::kmer_type Kmer_t;

  const size_t nh = kmap.size_hash();
  size_t n_aff = 0;
  size_t n_aff_exist = 0;
  for (size_t ih = 0; ih != nh; dots_pct(ih++, nh)) {

    const KMER_REC_t & krec0 = kmap[ih];
    if (krec0.is_valid_kmer()) { 
      unsigned n_pres = krec0.n_prefixes();
      unsigned n_sufs = krec0.n_suffixes();
      
      // ---- look at prefixes

      for (unsigned i = 0; i != n_pres; i++) {
	n_aff++;
	KmerFWRC<Kmer_t> kmerFR(krec0);
	kmerFR.push_left(krec0.prefix(i));
	if (kmap(kmerFR.canonical()).is_valid_kmer())
	  n_aff_exist++;
      }

      for (unsigned i = 0; i != n_sufs; i++) {
	n_aff++;
	KmerFWRC<Kmer_t> kmerFR(krec0);
	kmerFR.push_right(krec0.suffix(i));
	if (kmap(kmerFR.canonical()).is_valid_kmer())
	  n_aff_exist++;
      }
      
    }
  }      
    
  cout << "n_affixes= " << n_aff << endl
       << "n_found  = " << n_aff_exist << endl;
}





















template<class KMER_REC_t>
void base_vec_extension(const unsigned K, 
                        const KmerMap<KMER_REC_t> & kmap,
                        const BaseVec & bv,
                        const unsigned nb_extra,
                        BaseVec * bv_extra_p,
                        const bool get_pre)
{
  typedef typename KMER_REC_t::kmer_type Kmer_t;

  const size_t nb = bv.size();
  bv_extra_p->clear();

  SubKmers<BaseVec, Kmer_t> sub_kmer0(K, bv, get_pre ? 0 : nb - K);
  const Kmer_t kmer0 = get_pre ? sub_kmer0.rc() : sub_kmer0.fw();

  for (size_t i = 0; i < K; i++)
    bv_extra_p->push_back(kmer0[i]);
  
  KmerAffixesMapNavigator<KMER_REC_t> knav(kmap, kmer0);

  // get_pre = true  =>  follow suffixes of RC of first kmer
  // get_pre = false =>  follow suffixes of FW of last  kmer

  while (bv_extra_p->size() < K + nb_extra && 
         knav.n_suffixes() == 1) {
    bv_extra_p->push_back(knav.suffix(0));
    knav.next_by_suffix(0);
  }
}



template<class KMER_REC_t>
void base_vec_extensions_compute(const unsigned K, 
                                 const KmerMap<KMER_REC_t> & kmap,
                                 const BaseVecVec bvv_in,
                                 BaseVecVec * bvv_adj_p,
                                 const unsigned nb_extra)
{
  typedef typename KMER_REC_t::kmer_type Kmer_t;

  const size_t nbv = bvv_in.size();
  
  for (size_t ibv = 0; ibv < nbv; ibv++) {
    const BaseVec & bv_in = bvv_in[ibv];

    // ---- find prefix extension

    BaseVec bv_pre;
    base_vec_extension(K, kmap, bv_in, nb_extra, &bv_pre, true);

    if (bv_pre.size() == K + nb_extra)
      bvv_adj_p->push_back(bv_pre);


    // ---- find suffix extension

    BaseVec bv_suf;
    base_vec_extension(K, kmap, bv_in, nb_extra, &bv_suf, false);

    if (bv_suf.size() == K + nb_extra)
      bvv_adj_p->push_back(bv_suf);
  }
}












#endif
