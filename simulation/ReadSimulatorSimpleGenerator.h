///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

//
//  Author: Neil Weisenfeld - Sep 27, 2012 
//
// This is a class-ification of ReadSimulatorSimple.cc that works with Framework.h
//

#ifndef RSSG2_H_
#define RSSG2_H_

#include "simulation/Framework.h"
#include "simulation/ReadSimulatorSimpleCore.h"


////////////////////////////////////////////////////////////////////
// Top Level Generator: Reference+Reads+Errors
////////////////////////////////////////////////////////////////////


class ReadSimulatorSimpleGenerator : public Generator {
private:
  // (potential) output products
  BaseVecVec*		_bvv_p;
  QualVecVec*		_qvv_p;
  RefLocusVec*		_refv_p;
  ReadErrorVecVec*	_revv_p;
  PairsManager*		_pm_p;


  // parameters borrowed from the outside world
  const BaseVecVec& 	_ref;
  const vec<bool>	_ref_circular;

  unsigned 		_n_threads;

  void buildReadsFromPositions_parallel_helper(double err_del, double err_ins, double err_sub, const unsigned random_seed,
					size_t i_thread, size_t n_threads );

  size_t make_pair( size_t, size_t, size_t, size_t, size_t, bool);
  size_t make_read( size_t, size_t, size_t, bool );

  void compute_reads_size_and_position(bool paired, const size_t ref_id, const float coverage,
      const unsigned random_seed, const bool fw_only, const double prob_tile, RandomLength insert_len,

      const bool trimmable = false,  RandomLength read_len = RandomLength(0,0,0,0), const String& libname = "random_lib");

public:
  ReadSimulatorSimpleGenerator(
      const BaseVecVec& ref,
      const vec<bool>& ref_circular,
      const unsigned n_threads = 0,
      bool pairs = false,
      bool quals = false,
      bool readerrors=false,
      bool refloci = false) :
	_bvv_p(0),
	_qvv_p(0),
	_refv_p(0),
	_revv_p(0),
	_pm_p(0),
	_ref(ref),
	_ref_circular(ref_circular),
	_n_threads( boundNumThreads(n_threads) ) {

    /* always */	_bvv_p 	= new(BaseVecVec);
    if (quals) 		_qvv_p 	= new(QualVecVec);
    if (pairs) 		_pm_p 	= new(PairsManager);
    if (readerrors) 	_revv_p = new(ReadErrorVecVec);
    /* always */ 	_refv_p = new(RefLocusVec);
  };


  void computeReadsSizeAndPositionUnpaired(const size_t ref_id, const float coverage,
      const unsigned random_seed, const bool fw_only, const double prob_tile, RandomLength read_len,
      const bool trimmable = false);

  void computeReadsSizeAndPositionPaired(const size_t ref_id, const float coverage,
      const unsigned random_seed, const bool fw_only, const double prob_tile, RandomLength insert_len,
      RandomLength read_len, const bool trimmable = false, const String& libname = "random_lib");

  void buildReadsFromPositionsWithErrors(const unsigned random_seed, double err_del, double err_ins, double err_sub );

  void buildReadsFromPositions( const unsigned random_seed );	// separate interface because we may make this faster in the future

  ~ReadSimulatorSimpleGenerator() {
    delete _bvv_p;
    delete _qvv_p;
    delete _refv_p;
    delete _revv_p;
    delete _pm_p;
  };



  // output getters
  const BaseVecVec& getReads() const {
    ForceAssert(_bvv_p);
    return *_bvv_p;
  }

  const QualVecVec& getQuals() const {
    ForceAssert(_qvv_p);
    return *_qvv_p;
  }

  const RefLocusVec& getReadLocs() const {
    ForceAssert(_refv_p);
    return *_refv_p;
  }
  const ReadErrorVecVec& getReadErrors() const {
    ForceAssert(_revv_p);
    return *_revv_p;
  }

  const PairsManager& getPairsManager() const {
    ForceAssert(_pm_p);
    return *_pm_p;
  }
};


#endif /* RSSG2_H_ */
