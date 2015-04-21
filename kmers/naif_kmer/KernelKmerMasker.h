///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Aug 6, 2013 - <crdhelp@broadinstitute.org>
//


#ifndef KERNELKMERMASKER_H_
#define KERNELKMERMASKER_H_

#include "system/Assert.h"

// like a BV loc, but meant to be a genome location
class GVLoc
{
  uint64_t _ibv : 10;    // up to 2^10 = 1024 contigs
  uint64_t _ib  : 36;    // up to 2^36 base positions 68G bases
  uint64_t _dir :  2;    // 0:invalid  1:fw  2:rc  3:palindrome

public:
  GVLoc() : _ibv(0), _ib(0), _dir(0) {}

  void set_ibv       (const size_t ibv) { _ibv = ibv; }
  void set_ib        (const size_t ib)  { _ib  = ib; }
  void set_fw        (const bool fw)    { _dir = (fw ? 1 : 2); }
  void set_rc        (const bool rc)    { _dir = (rc ? 2 : 1); }
  void set_palindrome(const bool pal)   { if (pal) _dir = 3; }

  size_t ibv()           const { return _ibv; }
  size_t ib()            const { return _ib; }
  bool   is_fw()         const { return _dir == 1; }
  bool   is_rc()         const { return _dir == 2; }
  bool   is_palindrome() const { return _dir == 3; }
  bool   is_valid_loc()  const { return _dir != 0; }

  String to_str() const
  {
    return (String((is_palindrome() ? "dir= pal" :
                    is_fw() ? "dir= fw " :
		    is_rc() ? "dir= rc " : "dir= bad")) +
            " ibv/ib=" + ToString(_ibv) +
            " " + ToString(_ib));
  }

  friend
  bool operator < (const GVLoc & a, const GVLoc & b)
  {
    if (a._ibv < b._ibv) return true;
    if (a._ibv > b._ibv) return false;
    return (a._ib < b._ib);
  }
};




template<class KMER_t>
class KmerGVLoc : public KMER_t, public GVLoc
{
public:
  KmerGVLoc(const KMER_t & kmer) : KMER_t(kmer), GVLoc() {}
  explicit KmerGVLoc(const unsigned K = 0) : KMER_t(K), GVLoc() {}

  friend
  bool operator < (const KmerGVLoc & a, const KmerGVLoc & b)
  {
    if (static_cast<const KMER_t &>(a) < static_cast<const KMER_t &>(b)) return true;
    if (static_cast<const KMER_t &>(b) < static_cast<const KMER_t &>(a)) return false;
    return (static_cast<const GVLoc &>(a) < static_cast<const GVLoc &>(b));
  }

};





// ---- KernelKmerMasker -- build a mask of kmer positions passing the validator
//
template <class ELEM_t>
class KernelKmerMasker
{
private:
    const BaseVecVec   & _bvv;
    const size_t         _K;

    vecbitvector*        _p_bitvec;
    LockedData           _lock;     // lock for merge()

    vec<GVLoc>		 _baserange_tmp;

    const Validator*	 _p_validator;

    typedef typename ELEM_t::kmer_type  kmer_type;

public:
    typedef KmerGVLoc<kmer_type> rec_type;		// force KmerLcLoc records

    KernelKmerMasker(const BaseVecVec & bvv,
		     const size_t       K,
		     vecbitvector*      p_bitvec,
		     const Validator*   p_validator) :
	_bvv(bvv),
	_K(K),
	_p_bitvec(p_bitvec),
	_lock(),
	_baserange_tmp(),
	_p_validator(p_validator)
    {
	ForceAssert(p_bitvec);
	Mimic(_bvv, *_p_bitvec);
	for ( auto& bv : *_p_bitvec ) bv.Zero();
    };


    // copy constructor for temporary kernels
    explicit KernelKmerMasker(const KernelKmerMasker<ELEM_t> & that) :
	_bvv(that._bvv),
	_K(that._K),
	_p_bitvec(0),
	_lock(),
	_p_validator(that._p_validator)
	{}

    // interface function needed by naif_kmerize()
    size_t K() const { return _K; }

    // interface function needed by naif_kmerize()
    const BaseVecVec & bases() const { return _bvv; }


    // interface function needed by naif_kmerize()
    void parse_base_vec(ParcelBuffer<rec_type> * p_buf,
		      const size_t ibv)
    {
	ParcelBuffer<rec_type> & parcels = *p_buf;
	SubKmers<BaseVec, rec_type> kmer_cur(_K, _bvv[ibv]);
	while (kmer_cur.not_done()) {
	  rec_type & kmer = kmer_cur.canonical();
	  if (parcels.in_one_parcel(kmer)) {
	    kmer.set_ibv(ibv);
	    kmer.set_ib(kmer_cur.index_start());

	    if (kmer_cur.is_palindrome())  kmer.set_palindrome(true);
	    else                           kmer.set_fw(kmer_cur.is_canonical_fw());

	    parcels.add(kmer);
	  }
	  kmer_cur.next();
	}
    }


    // interface function needed by naif_kmerize()
    void summarize(const vec<rec_type> & krecs,
		   const size_t i0k,
		   const size_t i1k)
    {
	const size_t kf = i1k - i0k;

	GVLoc tmp;
	if ( !_p_validator || (*_p_validator)(kf) ) {
	    for ( size_t i = i0k; i < i1k; ++i ) {
		tmp.set_ib( krecs[i].ib() );
		tmp.set_ibv( krecs[i].ibv() );
		// the following is not needed, but here in case anyone
		// decides to hack the code for other purposes
		if ( krecs[i].is_palindrome() ) tmp.set_palindrome( true );
		else tmp.set_fw( krecs[i].is_fw() );
		_baserange_tmp.push_back(tmp);
	    }
	}
    }


    // interface function needed by naif_kmerize()
    void merge(const KernelKmerMasker<ELEM_t> & kernel_tmp,
	       const size_t i_parcel)
    {
	Locker lock(_lock);
	size_t K = kernel_tmp._K;
	for ( auto loc : kernel_tmp._baserange_tmp ) {
	    size_t ib  = loc.ib();
	    size_t ibv = loc.ibv();

	    ForceAssertLe(ibv, _p_bitvec->size() );
	    ForceAssertLe(ib + K - 1, (*_p_bitvec)[loc.ibv()].size() );

	    for ( size_t i = 0; i < K; ++i )
		(*_p_bitvec)[ibv].Set(ib, True);
	}
    }

};




#endif /* KERNELKMERMASKER_H_ */
