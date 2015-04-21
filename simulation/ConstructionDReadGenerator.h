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

#include "AlignmentProfile.h"
#include "Basevector.h"
#include "Bitvector.h"
#include "FeudalMimic.h"
#include "MemberOf.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "ReadLocationLG.h"
#include "TokenizeString.h"
#include "feudal/BinaryStream.h"
#include "kmers/KmerShape.h"
#include "math/Functions.h"
#include "paths/simulation/ErrorGenerator.h"
#include "random/NormalRandom.h"
#include "random/Random.h"
#include "random/Shuffle.h"
#include "solexa/SolexaMetrics.h" // solexa_metric_db
#include "system/Properties.h"

#include "simulation/Framework.h"
#include "simulation/ReadSimulatorSimpleCore.h"
#include "paths/simulation/ErrorGenerator.h"

////////////////////////////////////////////////////////////////////
// Top Level Generator: Reference+Reads+Errors
////////////////////////////////////////////////////////////////////


class ConstructionDReadGenerator : public Generator {
private:
  // (mandatory) output products
  BaseVecVec*		_bvv_p;		// we always generate the reads
  QualVecVec*		_qvv_p;		// we always generate the quals

  // (potential) output products -- see constructor (out_perfect, out_loci, out_pairs)
  BaseVecVec*		_perfect_p;	// perfect reads -- only made if not null
  RefLocusVec*		_refv_p;	// loci on the reference for each read -- only made if not null
  PairsManager*		_pm_p;		// pairs manager -- only made if not null; populated if we make pairs 

  // parameters borrowed from the outside world
  bool			_revcomp;	// whether to generate reads from rc inserts
  const ErrorGenerator* _errgen_p;	
  const BaseVecVec& 	_ref;

  size_t 		_n_threads;

  // internal per-thread helpers
  void CreateRandomReads( unsigned int nreads,  int n, size_t i_thread );
  void CreateRandomPairs( unsigned int npairs, int N, int n, const double dev = 0, const int jmean = -1, const double jdev = -1, size_t i_thread = 0);

  // this is the workhorse that builds pairs (when jmean,jdev=-1,-1) or jumps
  void buildPairs( size_t read_size, size_t coverage, size_t insert_size, double dev, int jmean, double jdev );

public:
  ConstructionDReadGenerator(
      const BaseVecVec& ref,
      bool revcomp = true,
      const ErrorGenerator* errgenp = 0,
      const size_t n_threads = 0,
      bool out_perfect = false,
      bool out_loci = false,
      bool out_pairs = false) :
	_bvv_p(0),
	_qvv_p(0),
	_perfect_p(0),
	_refv_p(0),
	_pm_p(0),
	_revcomp(revcomp),
	_errgen_p(errgenp),
	_ref(ref),
	_n_threads( boundNumThreads(n_threads) ) {

    ForceAssert(_errgen_p);

    /* always */		_bvv_p 	= new(BaseVecVec);
    /* always */ 		_qvv_p 	= new(QualVecVec);
    if (out_perfect)		_perfect_p = new(BaseVecVec);
    if (out_loci) 		_refv_p = new(RefLocusVec);
    if (out_pairs) 		_pm_p 	= new(PairsManager);
  };



  ~ConstructionDReadGenerator() {
    delete _bvv_p;
    delete _qvv_p;
    delete _pm_p;
    delete _perfect_p;
    delete _refv_p;
  };


  void buildUnpairedReads( size_t read_size, size_t coverage );
  void buildPairedReads( size_t read_size, size_t coverage, size_t insert_size, double dev );
  void buildPairedJumps( size_t read_size, size_t coverage, size_t insert_size, double dev, int jmean, double jdev );

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

  const PairsManager& getPairsManager() const {
    ForceAssert(_pm_p);
    return *_pm_p;
  }

  const BaseVecVec& getPerfect() {
    ForceAssert(_perfect_p);
    return *_perfect_p;
  }
};


#endif /* RSSG2_H_ */
