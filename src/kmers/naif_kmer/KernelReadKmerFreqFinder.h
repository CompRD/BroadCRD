///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#ifndef KMERS__NAIF_KMER__KERNEL_READ_KMER_FREQ_FINDER_H
#define KMERS__NAIF_KMER__KERNEL_READ_KMER_FREQ_FINDER_H

#include "kmers/naif_kmer/LockedBlocks.h"

#include "kmers/naif_kmer/Kmers.h"




struct KFStats 
{ 
  size_t max, min, sum; 
  KFStats() : max(0), min(0), sum(0) {}
};
  

template<class KMER_t>
class ReadKmerFreqFinder
{
  const size_t             _n_threads;
  const size_t             _n_blocks;
  LockedBlocks             _blocks;
  LockedBlocks           * _p_blocks;
  
  // local buffers
  //    first: i_read    
  //    second: kf 
  vec< vec< pair<size_t, size_t> > > _i_read_kf;  

  const size_t             _K;
  const BaseVecVec       & _bases;
 
  vec<KFStats>           & _kf_stats;

public:

  typedef  KmerBVLoc<KMER_t>   rec_type;

  ReadKmerFreqFinder(const size_t       K,
		     const BaseVecVec & bases,
		     vec<KFStats>     * p_kf_stats,
		     const size_t       n_threads)
    : _n_threads(n_threads),
      _n_blocks(2 * n_threads),
      _blocks(_n_blocks),
      _p_blocks(&_blocks),
      _i_read_kf(_n_blocks),
      _K(K),
      _bases(bases),
      _kf_stats(*p_kf_stats)
  {
    ForceAssertGt(_n_threads, 0ul);
  }
    

  ~ReadKmerFreqFinder()
  {}
      

  
  // copy constructor for temporary kernels
  explicit ReadKmerFreqFinder(const ReadKmerFreqFinder & that)
    : _n_threads(that._n_threads),
      _n_blocks(that._n_blocks),
      _blocks(),                   // no local blocks
      _p_blocks(that._p_blocks),   // get the blocks from 'that'
      _i_read_kf(_n_blocks),
      _K(that._K),
      _bases(that._bases),
      _kf_stats(that._kf_stats)
  {
    ForceAssertGt(_n_threads, 0ul);
  }
 
  // interface function needed by naif_kmerize()
  size_t K() const { return _K; }

  // interface function needed by naif_kmerize()
  const BaseVecVec & bases() const { return _bases; }
  


  // interface function needed by naif_kmerize()
  void parse_base_vec(ParcelBuffer<rec_type> * p_buf, 
                      const size_t ibv)
  {
    ParcelBuffer<rec_type> & parcel_buf = *p_buf;
    SubKmers<BaseVec, rec_type> kmer_cur(_K, _bases[ibv]);
    while (kmer_cur.not_done()) {
      rec_type & kmer = kmer_cur.canonical();
      if (parcel_buf.in_one_parcel(kmer)) {
        kmer.set_ibv(ibv);          // just need the i_read (or ibv) for kf stats
        parcel_buf.add(kmer);
      }
      kmer_cur.next();
    }
  } 



  void flush_bufs(const size_t max_buf_size)
  {
    bool verbose = false;//(max_buf_size == 0);

    std::set<size_t> blocks_todo;
    for (size_t i = 0; i < _n_blocks; i++)
      if (_i_read_kf[i].size() > max_buf_size)  // buffers too big? flush them!
        blocks_todo.insert(i);

    
    
    while (blocks_todo.size() > 0) {

      // ---- lock a subset of the data to update

      const size_t i_blk = _p_blocks->lock_some_block(blocks_todo);

      

      // ---- set global _kf_stats[i_read]

      vec< pair<size_t, size_t> > & pairs = _i_read_kf[i_blk];
      const size_t n = pairs.size();
      for (size_t i = 0; i < n; i++) {
	const size_t i_read = pairs[i].first;
	const size_t kf     = pairs[i].second;
	if (_kf_stats[i_read].sum == 0) {
	  _kf_stats[i_read].max = _kf_stats[i_read].min = kf;
	}
	else {
	  if (_kf_stats[i_read].max < kf) _kf_stats[i_read].max = kf;
	  if (_kf_stats[i_read].min > kf) _kf_stats[i_read].min = kf;
	}
	_kf_stats[i_read].sum += kf;
      }
      pairs.clear();



      // ---- unlock block and remove it from to-do list

      _p_blocks->unlock_block(i_blk);
      blocks_todo.erase(i_blk);

    } // while (blocks_todo.size() > 0) 

  } // void flush_bufs(...)



  size_t i_block(const size_t i, const size_t n) 
  { 
    return (i * _n_blocks) / n;  
  }


  // ---- summarizes all data in a tmp kernel.
  // interface function needed by naif_kmerize()
  void summarize(const vec<rec_type> & parcel, 
                 const size_t i0k,
                 const size_t i1k)
  {
    const size_t kf = i1k - i0k;
    const size_t n_reads = _bases.size();
    for (size_t ik = i0k; ik != i1k; ik++) {
      const size_t i_read = parcel[ik].ibv();
      const size_t i_blk  = i_block(i_read, n_reads);
      _i_read_kf[i_blk].push_back(pair<size_t, size_t>(i_read, kf));
    }
    flush_bufs(5000);
  }
  
 










  // ---- this merges all the buffered data in tmp kernels.
  // interface function needed by naif_kmerize()
  void merge(ReadKmerFreqFinder & tmp_kernel,
	     const size_t i_parcel)
  {
    // final flush 
    tmp_kernel.flush_bufs(0);
  }






};  // class ReadKmerFreqFinder
  







#endif
