///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

//
//  Author: Neil Weisenfeld - Oct 11, 2012
//

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "simulation/ReadSimulatorSimpleGenerator.h"
#include "simulation/ReferenceIterator.h"

size_t ReadSimulatorSimpleGenerator::make_pair(size_t nb_insert, size_t len1, size_t len2, size_t ref_id, size_t ib_ref, bool rc)
{

  RefLocus read1( ref_id, ib_ref, len1, false );
  RefLocus read2( ref_id, ib_ref + nb_insert - len2, len2, true);

  if ( !rc ) {
    _refv_p->push_back( read1 );
    _refv_p->push_back( read2 );

  } else {
    _refv_p->push_back( read2 );
    _refv_p->push_back( read1 );
  }

  return len1+len2;

}


size_t ReadSimulatorSimpleGenerator::make_read(size_t nb_read, size_t ref_id, size_t ib_ref, bool rc )
{
  _refv_p->push_back(
 	      RefLocus( ref_id, ib_ref, nb_read, rc)
 	      );

  return nb_read;
}


void ReadSimulatorSimpleGenerator::compute_reads_size_and_position(bool paired, const size_t ref_id, const float coverage,
      const unsigned random_seed, const bool fw_only, const double prob_tile,
      RandomLength insert_len, const bool trimmable, RandomLength read_len, const String& libname)
{
  const size_t nb_ref = _ref[ref_id].size();
  const size_t nb_reads_max = nb_ref * coverage;

  const bool tiled = ( prob_tile > 0.0 );

  RandomGen random(random_seed + 2);

  size_t nb_reads = 0;


  size_t jb_ref = 0;
  size_t ib0_ref = 0;   // the starting index in the reference base vec

  const bool is_bv_circular = _ref_circular[ref_id];

  while (nb_reads < nb_reads_max) {

    if (nb_reads == 0 ||     // 1st time
        jb_ref >= nb_ref) {  // reached the end
      jb_ref = 0;
      ib0_ref = (is_bv_circular ? random.float01() * nb_ref : 0);
    }


    size_t nb_read = insert_len.value();  // get a random read length
    size_t tile_skip = nb_read;		// by default we skip a full read if we're tiling

    // prob_tile == 1.0 => tiled reads
    // prob_tile == 0.5 => only half of tiled reads
    if (!tiled || random.float01() < prob_tile) { // include read

      size_t ib_ref = 0;

      if (tiled) {
	ib_ref = (ib0_ref + jb_ref) % nb_ref;

	if (jb_ref + nb_read  > nb_ref) { // trim read length to fit?
	 // needs trimming
	  if ( trimmable )
	    nb_read = nb_ref - jb_ref;
	  else
	    nb_read = 0;	// if we can't trim, but need to, we'll skip it
	}
      }
      else { // not tiled
        ib_ref = (is_bv_circular ? random.unsignedN(nb_ref) : random.unsignedN(nb_ref - nb_read));
        //cout << "ib_ref= " << ib_ref << endl;
      }

      if ( nb_read > 0 ) {

	if ( paired ) {
	  size_t len1 = read_len.value();
	  size_t len2 = read_len.value();

	  ForceAssertLe(len1, nb_read);
	  ForceAssertLe(len2, nb_read);

	  nb_reads += make_pair( nb_read, len1, len2, ref_id, ib_ref,(fw_only ? false : random.unsignedN(2)));

	  if ( _pm_p ) {
	    size_t id2 = _refv_p->size();
	    int sep = insert_len.mu() - 2 * read_len.mu();
	    _pm_p->addPair(id2-2, id2-1, sep, insert_len.sig(), libname, true );
	  }

	  tile_skip = len1;

	} else {

	  nb_reads += make_read( nb_read, ref_id, ib_ref, (fw_only ? false : random.unsignedN(2)));
	  tile_skip = nb_read;

	}

		// for coverage purposes, how many bases did we cover -- don't include padding

      } // if nb_read > 0
    } // if tiled read accepted

    jb_ref += tile_skip;
  }
}

void ReadSimulatorSimpleGenerator::computeReadsSizeAndPositionUnpaired(const size_t ref_id, const float coverage,
      const unsigned random_seed, const bool fw_only, const double prob_tile,
      RandomLength read_len, const bool trimmable) {
  compute_reads_size_and_position(false, ref_id, coverage, random_seed, fw_only, prob_tile, read_len, trimmable);
}

void ReadSimulatorSimpleGenerator::computeReadsSizeAndPositionPaired(const size_t ref_id, const float coverage,
      const unsigned random_seed, const bool fw_only, const double prob_tile, RandomLength insert_len,
      RandomLength read_len, const bool trimmable, const String& libname ) {
  compute_reads_size_and_position(true, ref_id, coverage, random_seed, fw_only, prob_tile, insert_len, trimmable,read_len, libname);
}

void ReadSimulatorSimpleGenerator::buildReadsFromPositions_parallel_helper( double err_del, double err_ins, double err_sub, const unsigned random_seed, const size_t i_thread, const size_t n_threads )
{
  RandomGen random(random_seed + 3 + i_thread); // a generator for this thread, no need to lock. faster this way!

//  const size_t nbv_ref = bvv_ref.size();
  size_t nbv_tot = _bvv_p->size();
  size_t ibv0, ibv1, nbv;

  AssertGt( nbv_tot, n_threads );	// just to be sure

  ibv0 =  i_thread     * nbv_tot / n_threads;
  ibv1 = (i_thread + 1)* nbv_tot / n_threads;

  nbv = ibv1 - ibv0;

  for (size_t jbv = 0; jbv < nbv; jbv++) {
    // if (i_thread == 0) dots_pct(jbv, nbv);

    const size_t ibv = ibv0 + jbv;

    RefLocus& rl = (*_refv_p)[ibv];

    const BaseVec & bv_ref = _ref[rl.mRefID];
    const size_t nb_ref = bv_ref.size();

    const size_t nb = rl.mLength;
    const size_t ib0 = rl.mOffset + ( rl.mRC ? rl.mLength - 1 : 0 );

    BaseVec & bv = (*_bvv_p)[ibv];
    bv.resize(nb);

    if ( _qvv_p )
      (*_qvv_p)[ibv].resize(nb);

    ReferenceIterator ref_iter( bv_ref, ib0, rl.mRC, _ref_circular[rl.mRefID], err_del, err_ins, err_sub, &random );
    size_t ib = 0;
    for ( ib = 0; ib < nb && !ref_iter.done() ; ib++, ++ref_iter )
      bv.set(ib, *ref_iter );

    if ( ib < nb ) {
      std::cout << "ib=" << ib << ", nb=" << nb << ", done=" << ref_iter.done() << ", err_sum=" << err_del+err_ins+err_sub << std::endl;
    }

    const ReadErrorVec rev = ref_iter.getErrors();

    if ( _qvv_p ) {
      QualVec& qv = (*_qvv_p)[ibv];
      for ( size_t i = 0; i < nb; ++i ) {		// model increasing drop with position from Q40 to Q20
	qv[i] = log(static_cast<float>(nb - i)+1.0) / log( nb + 1.0) * (40U-20U) + 20U;
      }
      for ( size_t i = 0; i < rev.size(); ++i ) {		// drop by 60% at the site of errors
	size_t loc = rev[i].getLocation();
	qv[loc] *= 0.6;
      }

    }

    if ( _revv_p )
      (*_revv_p)[ibv] = rev;

  }
}



void ReadSimulatorSimpleGenerator::buildReadsFromPositions(const unsigned random_seed)
{
  this->buildReadsFromPositionsWithErrors(random_seed,0.,0.,0.);
}


void ReadSimulatorSimpleGenerator::buildReadsFromPositionsWithErrors(const unsigned random_seed, double err_del, double err_ins, double err_sub )
{
  // initialize vectors which will have members modified in parallel
  const size_t n_reads = _refv_p->size();
  ForceAssertGt( n_reads, 0U );	// must have already made the reads
  _bvv_p->resize( n_reads );

  if ( _revv_p )
  _revv_p->resize(n_reads);

  if ( _qvv_p )
    _qvv_p->resize(n_reads);

#pragma omp parallel for
  for (size_t i_thread = 0; i_thread < _n_threads; i_thread++ )
    buildReadsFromPositions_parallel_helper(err_del,err_ins,err_sub,random_seed,i_thread,_n_threads);

}

