///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////



#ifndef KMERS__NAIF_KMER__UNIBASES__KMERS_H
#define KMERS__NAIF_KMER__UNIBASES__KMERS_H





class UnibaseNode
{
public:
  vec<size_t> iub_pre;
  vec<size_t> iub_suf;
  vec<bool> fw_pre;
  vec<bool> fw_suf;

  template<class KMER_REC_t>
  UnibaseNode(const size_t K, 
              const BaseVec & bv_ub,
              const KmerMap<KMER_REC_t> & kmap)
  {
    _find_prefixes(K, bv_ub, kmap);
    _find_suffixes(K, bv_ub, kmap);
  }

private:
  template<class KMER_REC_t>
  void _find_prefixes(const size_t K, 
		      const BaseVec & bv_ub,
		      const KmerMap<KMER_REC_t> & kmap)
  {
    _find_affixes(K, bv_ub, kmap, true,  & iub_pre, & fw_pre);
  }


  template<class KMER_REC_t>
  void _find_suffixes(const size_t K, 
		      const BaseVec & bv_ub,
		      const KmerMap<KMER_REC_t> & kmap)
  {
    _find_affixes(K, bv_ub, kmap, false, & iub_suf, & fw_suf);
  }


  template<class KMER_REC_t>
  void _find_affixes(const size_t K, 
		     const BaseVec & bv_ub,
		     const KmerMap<KMER_REC_t> & kmap,
		     const bool get_pre,
		     vec<size_t> * iub_p,
		     vec<bool> * fw_p)
  {
    typedef typename KMER_REC_t::kmer_type Kmer_t;

    const size_t nb = bv_ub.size();
    
    SubKmers<BaseVec, Kmer_t> kmer0(K, bv_ub, get_pre ? 0 : nb - K);
    const bool       is_fw0 = kmer0.is_canonical_fw();
    
    const KMER_REC_t & krec0 = kmap(kmer0.canonical());
    const unsigned     n_aff = (get_pre ? 
                                krec0.n_prefixes(is_fw0) :
                                krec0.n_suffixes(is_fw0));
    
    // ---- go through all affixes
    
    for (unsigned i_aff = 0; i_aff < n_aff; i_aff++) {
      
      // ---- build the kmer + affix
      
      KmerFWRC<Kmer_t> kmer(kmer0.fw());
      if (get_pre) kmer.push_left (krec0.prefix(i_aff, is_fw0));
      else         kmer.push_right(krec0.suffix(i_aff, is_fw0));
      
      const KMER_REC_t & krec = kmap(kmer.canonical());
      const bool is_fw_ub     = krec.is_fw();           // if the unibase kmer is canonical
      const bool is_fw        = kmer.is_canonical_fw(); // if the affix  kmer is canonical
      // ---- save unibase index and whether they face the same way
      
      iub_p->push_back(krec.ibv()); 
      fw_p->push_back(is_fw == is_fw_ub);
    }
  }

};



template<class KMER_REC_t>
class UnibaseGraph
{
  const size_t                 _K;
  const KmerMap<KMER_REC_t>  & _kmap;
  const BaseVecVec           & _bvv_ub;
  const vec<float>           & _mean_kf;
  vec<UnibaseNode>             _graph;



public:
  UnibaseGraph(const size_t K,
               const KmerMap<KMER_REC_t> & kmap,
               const BaseVecVec & bvv_ub,
               const vec<float> & mean_kf) :
    _K(K),
    _kmap(kmap),
    _bvv_ub(bvv_ub),
    _mean_kf(mean_kf),
    _graph()
  {
    const size_t nub = _bvv_ub.size();
    ForceAssertEq(nub, _mean_kf.size());

    for (size_t iub = 0; iub < nub; iub++)
      _graph.push_back(UnibaseNode(K, _bvv_ub[iub], kmap));
  }

  unsigned n_prefixes(const size_t iub, const bool fw) const
  { return (fw ? _graph[iub].iub_pre.size() : _graph[iub].iub_suf.size()); }

  unsigned n_suffixes(const size_t iub, const bool fw) const 
  { return (fw ? _graph[iub].iub_suf.size() : _graph[iub].iub_pre.size()); }

  size_t iub_prefix(const size_t iub, const size_t i_pre, const bool fw) const 
  { return (fw ? _graph[iub].iub_pre[i_pre] : _graph[iub].iub_suf[i_pre]); }

  size_t iub_suffix(const size_t iub, const size_t i_suf, const bool fw) const
  { return (fw ? _graph[iub].iub_suf[i_suf] : _graph[iub].iub_pre[i_suf]); }

  bool fw_prefix(const size_t iub, const size_t i_pre, const bool fw) const 
  { return (fw ? _graph[iub].fw_pre[i_pre] : _graph[iub].fw_suf[i_pre]); }

  bool fw_suffix(const size_t iub, const size_t i_suf, const bool fw) const 
  { return (fw ? _graph[iub].fw_suf[i_suf] : _graph[iub].fw_pre[i_suf]); }

  float mean_kf(const size_t iub) const
  { return _mean_kf[iub]; }
    
  BaseVec & bv(const size_t iub) const 
  { return _bvv_ub[iub]; }

  size_t nk(const size_t iub) const 
  { return _bvv_ub[iub].size() - _K + 1; }


};
















template<class KMER_REC_t>
void kmer_extend_forward(const typename KMER_REC_t::kmer_type & kmer0,
			 KmerMap<KmerBVLoc<KMER_REC_t> > * kmap_p, // updated kmer locs
                         BaseVecVec * bvv_ub_p,
                         vec<float> * mean_kf_p = 0)
{
  typedef typename KMER_REC_t::kmer_type Kmer_t;
  const unsigned K = kmer0.K();
  const size_t ibv = bvv_ub_p->size();

  size_t nk = 0;
  size_t freq_sum = 0;
  unsigned base_suf = 0;
  BaseVec bv_ub;

  // ---- declare the walking kmer
  KmerFWRC<Kmer_t> kmerFR(kmer0);
  bool done = false;
  while (!done) {
    
    KmerBVLoc<KMER_REC_t> & krec = (*kmap_p)(kmerFR.canonical());
    if (! krec.is_valid_kmer()) {
      cout << "unexpected invalid kmer:" << endl;
      cout << "fw: " << hieroglyphs(kmerFR.fw()) << endl;
      cout << "rc: " << hieroglyphs(kmerFR.rc()) << endl;
      cout << "ibv= " << ibv << endl;
      cout << "nk= " << nk << endl;
      cout << "bv: " << hieroglyphs(bv_ub) << endl;
      abort();
    }
    done = true;

    if (! krec.tag()) {  // kmer hasn't bee visited
      const bool is_fw = kmerFR.is_canonical_fw();
      const bool is_pal = kmerFR.is_palindrome();

      if (nk == 0 ||                      // 1st kmer is always included
	  krec.n_prefixes(is_fw) == 1) {  // but after that only 1 prefix allowed
	
	krec.tag_set();      // kmer was visited
	krec.set_ibv(ibv);
	krec.set_ib(nk);
	if (is_pal) krec.set_palindrome(true);
	else        krec.set_fw(is_fw);

	// ---- update basevec
	if (nk == 0) 
	  for (unsigned i = 0; i < K; i++) 
	    bv_ub.push_back(kmer0[i]); 
	else 
	  bv_ub.push_back(base_suf);

	nk++;
	freq_sum += krec.freq(); 

	if (krec.n_suffixes(is_fw) == 1) {
	  base_suf = krec.suffix(0, is_fw);
	  kmerFR.push_right(base_suf);
	  done = false;
	}
      }
    }
  }
   
  if (bv_ub.size() > 0) {
    bvv_ub_p->push_back(bv_ub);
    if (mean_kf_p != 0)
      mean_kf_p->push_back(float(freq_sum) / float(nk));
  }
}


template<class KMER_REC_t>
void kmer_extend_backward(const typename KMER_REC_t::kmer_type kmer,
			  KmerMap<KMER_REC_t> * kmap_p, // updated with kmer locs
                          BaseVecVec * bvv_ub_p,
                          vec<float> * mean_kfs_p = 0)
{
  kmer_extend_forward(reverse_complement(kmer), kmap_p, bvv_ub_p, mean_kfs_p);
}










template<class KMER_REC_t>
void unibases_from_kmer_affix_map(const unsigned K, 
                                  KmerMap<KMER_REC_t> * kmap_p,
                                  BaseVecVec * bvv_ub_p,
                                  vec<float> * mean_kf_p = 0)
{
  typedef typename KMER_REC_t::kmer_type Kmer_t;

  const size_t nh = kmap_p->size_hash();

  // ---- Find all non 1-1 kmers and follow them

  cout << "Finding unibases:" << endl;

  for (size_t ih = 0; ih < nh; dots_pct(ih++, nh)) {

    const KMER_REC_t & krec = (*kmap_p)[ih];
    if (krec.is_valid_kmer()) {
      const size_t n_pres = krec.n_prefixes();
      const size_t n_sufs = krec.n_suffixes();

      if (n_sufs != 1 || n_pres != 1) { 
	const Kmer_t kmer = krec;
   
	// ---- follow self kmer
	
	if (n_pres != 1) { // 0, 2, 3, or 4
	  kmer_extend_forward(kmer, kmap_p, bvv_ub_p, mean_kf_p);
	}
	else if (n_sufs != 1) { // 0, 2, 3, or 4
	  kmer_extend_backward(kmer, kmap_p, bvv_ub_p, mean_kf_p);
	}

        // ---- follow the suffixes

	if (n_sufs > 1) {
	  for (unsigned i_suf = 0; i_suf < n_sufs; i_suf++) {
	    Kmer_t kmer_next = kmer;
	    kmer_next.push_right(krec.suffix(i_suf));
	    kmer_extend_forward(kmer_next, kmap_p, bvv_ub_p, mean_kf_p);
	  }
	}

        // ---- follow the prefixes

	if (n_pres > 1) {
	  for (unsigned i_pre = 0; i_pre < n_pres; i_pre++) {
	    Kmer_t kmer_next = kmer;
	    kmer_next.push_left(krec.prefix(i_pre));
	    kmer_extend_backward(kmer_next, kmap_p, bvv_ub_p, mean_kf_p);
	  }
	}

      }
    }

  }
  size_t nub = bvv_ub_p->size();
  
  cout << setw(10) << bvv_ub_p->size() << " unibases found." << endl;


  // ---- find all the 1-1 kmers left behind; jubt cycles

  cout << "Finding unibase cycles:" << endl;

  for (size_t ih = 0; ih < nh; dots_pct(ih++, nh)) {

    const KMER_REC_t & krec = (*kmap_p)[ih];
    if (krec.is_valid_kmer() && !krec.tag()) {

      const size_t n_pres = krec.n_prefixes();
      const size_t n_sufs = krec.n_suffixes();
      ForceAssertEq(n_pres, 1u);
      ForceAssertEq(n_sufs, 1u);
   
      // ---- follow 1-1 kmer; stops after 1 cycle period
      kmer_extend_forward(krec, kmap_p, bvv_ub_p, mean_kf_p);
    }

  }
  cout << setw(10) << bvv_ub_p->size() - nub << " unibase cycles found." << endl;
  cout << setw(10) << bvv_ub_p->size() << " unibases total found." << endl;

  

  // ---- Verify


  if (false) {
  
    size_t n_11 = 0;
    size_t n_11_bad = 0;
    size_t n_xx = 0;
    size_t n_xx_valid = 0;
    
    for (size_t ih = 0; ih < nh; dots_pct(ih++, nh)) {
      const KMER_REC_t & krec = (*kmap_p)[ih];
      if (krec.is_valid_kmer()) {
	
	const size_t n_pres = krec.n_prefixes();
	const size_t n_sufs = krec.n_suffixes();
	
	if (n_sufs == 1 && n_pres == 1) {
	  if (krec.tag())
	    n_11 ++;
	  else {
	    n_11_bad++;
	    //cout << hieroglyphs(krec) << endl;
	  }
	}
	else {
	  if (krec.tag())
	    n_xx_valid ++;
	  else 
	    n_xx++;
	}
      }
    }
    
    cout << "n_11=       " << n_11 << endl;
    cout << "n_11_bad=   " << n_11_bad << endl;
    cout << "n_xx=       " << n_xx << endl;
    cout << "n_xx_valid= " << n_xx_valid << endl;
  }
}
  






















#endif
