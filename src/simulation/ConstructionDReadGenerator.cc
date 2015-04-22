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

#include "simulation/ConstructionDReadGenerator.h"
#include "simulation/ReferenceIterator.h"
#include "RefLocus.h"


namespace {

/**
   Local function: CalculateGenomeBases

   Give a genome and a read or insert size, this function will calculate
   the number of distinct positions in the genome a simulated insert or read
   can be derived from. In other words, the number of possible start positions
   in the genome for a read or insert of the specified size.

   Input parameters:

     genome - the genome to use
     nbases - size of insert or read (in terms of bases on the reference)

   Output parameters:
     returns the number of genome bases as defined above

*/
longlong CalculateGenomeBases(const BaseVecVec& genome, const int nbases) {
  longlong gbases = 0;
  for ( size_t i = 0; i < genome.size( ); i++ ) {
    const int gsize = genome[i].isize();
    if ( gsize < nbases ) continue;
    gbases += gsize - nbases + 1;
  }
  return gbases;
}

/**
   Local function: FindContigAndUpdatePosition

   Related to <CalculateGenomeBases>, this function will find the contig and
   position on the contig, given a position on the genome. The position in the
   genome is defined in terms of the number of distinct positions in the genome
   a simulated insert or read can be derived from - as calculated in
   <CalculateGenomeBases>.

   Input parameters:

     genome - the genome to use
     nbases - size of insert or read (in terms of bases on the reference)
     pos - a position on the geneome

   Output parameters:
     returns the contig index where pos lies
     pos - the corresponding position on the contig

*/
genome_part_id_t FindContigAndUpdatePosition(const BaseVecVec& genome, const int nbases, longlong& pos) {
  for ( size_t i = 0; i < genome.size( ); i++ ) {
    const int gsize = genome[i].isize();
    if ( gsize < nbases ) continue;
    if ( pos > gsize - nbases )
      pos -= ( gsize - nbases + 1 );
    else {
      return i;
    }
  }
  FatalErr("Read start position in genome not valid");
  return 0; 	// notused -- this is to quell warnings from Eclipse about a missing return value.
}



long big_local_random(RNGen& rand) { return (rand.next()<<31) | rand.next();  }


} // namespace {}



/**
   Local function: CreateRandomReads

   Create simulated unpaired reads from a genome.

   Input parameters:

      genome - the known genome from which the simulated reads are generated
      n - read length; all simulated reads will have this exact length.
      nreads - the number of reads to generate; this has been calculated from
         genome size and the requested <coverage>.
      method - method to use to simulate reads with errors (default NONE)

   Output parameters:

      reads - the simulated reads.  note that these simulated reads will have
         artificially generated errors; true_reads gives the actual, error-free
         original simulated reads.
      true_reads - the true original simulated reads, before any errors are
         introduced.  only generated if `generate_true_reads' is true.
         of course, only for simulated reads can we know this with certainty.
      generate_true_reads - whether to generate true reads
      quals - qualities of the simulated reads
      pairs - the pairings of the simulated reads, showing for each simulated
         read pair which two reads in _reads_ comprise that pair
      loc_on_ref - the true location of each read on the genome (of course,
         only for simulated reads can we know this with certainty.)

   Shares code with <CreateRandomPairs()>.
*/
void
ConstructionDReadGenerator::
CreateRandomReads( unsigned int nreads,  int n, size_t i_thread)
{
  // thread housekeeping
  ForceAssertGe(_bvv_p->size(), nreads);
  size_t array_base = _bvv_p->size() - nreads;

  const size_t real_threads = std::min( static_cast<size_t>(nreads), _n_threads );

  if ( real_threads <= i_thread )
    return;	// no work to do

  RNGen local_random(i_thread);

  cout << "Creating " << nreads << " unpaired reads of length " << n << " ";


  // Determine number of bases in genome (gbases) a read can start from given
  // the minimum and maximum size (on the reference) of a read
  int min_bases = n, max_bases = n;
  if (_errgen_p)
    _errgen_p->GetMinMaxRefBaseCount(min_bases, max_bases, n);
  map<int, longlong> gbases;
  for ( int nbases = min_bases; nbases <= max_bases; ++nbases)
    gbases[nbases] = CalculateGenomeBases(_ref, nbases);

  // Generate each read in turn
  for ( size_t read = 0; read < nreads; read++ ) {
    size_t i = array_base + read;

    longlong pos = big_local_random( local_random );  // Random number used to pick pos on genome
    bool rc = local_random.next() & 1;

    int nbases = n; // number of bases on reference required to generate read

    AlignmentProfile profile;
    // Use error generator to pick error profile with quality scores
    if ( _errgen_p ) {
      qualvector qual;
      _errgen_p->GetErrorProfile(profile, qual, n);
      (*_qvv_p)[i] = qual;
      nbases = profile.GetRefBaseCount();  // number of bases on reference required
    }

    // Pick random position in entire genome for start of read
    pos %= gbases[nbases];

    // Find contig and position on contig for start of read
    int contig = FindContigAndUpdatePosition(_ref, nbases, pos);

    const basevector& g = _ref[contig];

    basevector b;
    b.SetToSubOf( g, pos, n );
    if (rc) b.ReverseComplement();
    if (_perfect_p) (*_perfect_p)[i]=b;

    // Use profile from error generator to create simulated reads with errors
    if ( _errgen_p ) {
      if ( rc ) {
	profile.CreateRcRead(g, pos+nbases-1, b);
      } else {
	profile.CreateRead(g, pos, b);
      }
    }

    if ( _refv_p )
      (*_refv_p)[i]=RefLocus(contig, pos, nbases, rc) ;
    (*_bvv_p)[i]=b;

  }
  cout << endl;
}



/**
   Local function: CreateRandomPairs

   Generate simulated paired reads from a genome. If npairs is large relative
   to the genome size, this will be slow.

   Input parameters:

      genome - the known genome from which the simulated reads are generated
      npairs - number of simulated read pairs to generate
      n - read length; all simulated reads will have this exact length.
      N - fragment length
      dev - the stddev for the variation in fragment size N (in the absolute
            number of bases, not as percentage of N); this simulates the fact
	    that in real <libraries>, there is variation in fragment sizes
            (although the laboratory procedures are tuned to try to minimize
	    this variation).
      jmean, jdev - the mean and stddev of the mini-fragment size in a jumping
                    library; otherwise -1 in which case this is not jumping
      method - method to use to simulate reads with errors (default NONE)

   Output parameters:

      reads - the simulated reads.  note that these simulated reads will have
         artificially generated errors; true_reads gives the actual, error-free
         original simulated reads.
      true_reads - the true original simulated reads, before any errors are
         introduced.  only generated if `generate_true_reads' is true.  of
	 course, only for simulated reads can we know this with certainty.
      generate_true_reads - whether to generate true reads
      quals - qualities of the simulated reads
      pairs - the pairings of the simulated reads, showing for each simulated
         read pair which two reads in _reads_ comprise that pair
      loc_on_ref - the true location of each read on the genome (of course,
         only for simulated reads can we know this with certainty.)

   Shares code with <CreateRandomReads()>.
*/
void ConstructionDReadGenerator::
CreateRandomPairs( unsigned int npairs, int N, int n,
		   const double dev, const int jmean,
		   const double jdev, size_t i_thread)
{
  /*
   * some thread housekeeping
   */
  ForceAssertGe(_bvv_p->size(), 2 * npairs);
  size_t array_base = _bvv_p->size() - ( 2 * npairs ); 

  const size_t real_threads = std::min( static_cast<size_t>(npairs), _n_threads );

  if ( real_threads <= i_thread )
    return;	// no work to do

  RNGen local_random(i_thread);

  // Normal distributions for insert size and jumping fragment size
  NormalRandom insert_normal( 0, dev );
  NormalRandom jfrag_normal( jmean, jdev );
  insert_normal.seed(i_thread);
  jfrag_normal.seed(i_thread);

  /* If this is a jumping library: Create a simulated jumping fragment.
   *
   * In real life, a jumping read pair is created as follows:  A long insert
   * (say, ~4 kbp) is tied together at its ends to create a loop; the tying
   * location is called the "cut point" and is marked with a biotin handle.
   * The loop is then sheared into small fragments.  The fragment containing
   * the handle is then read at both ends, like a normal insert, to create
   * the read pair.  Because the ends of this fragment come from opposite
   * ends of the original insert, the resulting reads are about 4kbp away
   * from each other, and are oriented in opposite directions.
   *
   * This can go wrong in 3 different ways, leading to 3 types of chimerism:
   * Read chimerism - One of the reads bridges the cut point.  That read
   *     will align correctly but partially to multiple places in the genome.
   * Fragment chimerism - The fragment that gets sequenced is an irrelevant
   *     fragment that does not contain the cut point, resulting in a pair
   *     that is perfectly normal except that it is non-jumping.
   * Insert chimerism - Instead of joining the ends of an insert, the biotin
   *     accidentally connects two inserts from different regions of the
   *     genome, resulting in a read pair with arbitrary (usually very large)
   *     distance and orientation.  May co-exist with read chimerism.
   *
   * The frequency (per pair) of read chimerism is 2n / jmean.
   * The frequencies of fragment and insert chimerism are fc_freq and ic_freq
   */

  bool jumping = ( jmean != -1 && jdev != -1 );

  // Generate each read pair in turn
  basevector frag;
  for ( size_t pair = i_thread; pair < npairs; pair+=real_threads ) {
    size_t i = array_base + ( 2*pair );

    // Insert length must be at least the size of a read
    int insert_length = N;
    if ( dev != 0 ) {
      insert_length += int( round( insert_normal.value( ) ) );
      insert_length = Max(insert_length, n);
    }
    // Determine number of bases in genome (posrange) an insert can start from,
    // given the insert size (including deviation)
    longlong posrange = CalculateGenomeBases(_ref, insert_length);

    // This could happen if insert_length is too big.
    if ( posrange < 1 )
      FatalErr( "\nFATAL ERROR!\nThe reference genome is too small for this library.\n" );

    // Pick random position in entire genome for start of insert
    longlong pos = big_local_random( local_random ) % posrange;

    // Find contig and position on contig for start of insert
    genome_part_id_t contig = FindContigAndUpdatePosition(_ref, insert_length, pos);

    const basevector& g = _ref[contig];
    if ( jumping ) {

      // Determine fragment length, based on the inputs jmean and jdev
      int frag_length = int( round( jfrag_normal.value( ) ) );
      frag_length = Min(frag_length, insert_length);
      // Fragment length must be at least 10 bases larger than read length.
      // This is because extra bases may need to be read off the fragment,
      // in case an error profile is used and it gives deletions.
      frag_length = Max(frag_length, n + 10);

      // Randomly determine fragment start/stop location on insert, and create
      // the fragment's basevector
      int frag_start, frag_end;
      basevector fragA, fragB;


      frag_end = big_local_random( local_random ) % frag_length;


      frag_start = frag_end - frag_length + insert_length;

      fragA.SetToSubOf( g, pos + frag_start, insert_length - frag_start );
      fragB.SetToSubOf( g, pos, frag_end );
      frag = Cat( fragA, fragB );

      // Set read locations

      if ( _revcomp == True && (local_random.next() & 0x01) ) {			// create approx 50% fragments from RC strand
	frag.ReverseComplement();
	if ( _refv_p ) {
	  (*_refv_p)[i] = RefLocus(contig, pos + frag_end - n, n, ORIENT_FW );
	  (*_refv_p)[i+1] = RefLocus(contig, pos + frag_start, n, ORIENT_RC );
	}
      } else if ( _refv_p ) {
	  (*_refv_p)[i] = RefLocus(contig, pos + frag_start, n, ORIENT_RC );
	  (*_refv_p)[i+1] = RefLocus(contig, pos + frag_end - n, n, ORIENT_FW );
      }
    }


    // If not a jumping library: Create a regular insert fragment
    else {
      frag.SetToSubOf( g, pos, insert_length );

      if (_revcomp == True && (local_random.next() & 0x01) ) {			// create approx 50% fragments from RC strand
	frag.ReverseComplement();
	if ( _refv_p ) {
	  (*_refv_p)[i] = RefLocus( contig, pos + insert_length - n, n, ORIENT_RC );
	  (*_refv_p)[i+1] = RefLocus( contig, pos, n, ORIENT_FW );
	}
      } else if ( _refv_p ) {
	(*_refv_p)[i] = RefLocus( contig, pos, n, ORIENT_FW );
	(*_refv_p)[i+1] = RefLocus( contig, pos + insert_length - n, n, ORIENT_RC );
      }
    }


    // Create the reads by reading this fragment
    int frag_length = frag.size( );
    basevector b1, b2;
    b1.SetToSubOf( frag, 0, n );
    b2.SetToSubOf( frag, frag_length - n, n );
    b2.ReverseComplement( );

    // If this is a jumping fragment, reverse the reads so that they will be
    // pointing toward each other on the reference
    // This step is performed manually on non-simulated jumping libraries and
    // is required for the insert-walking functions in the RunAllPaths pipeline
    // to work properly
    if ( jumping ) {
      b1.ReverseComplement( );
      b2.ReverseComplement( );
    }


    if (_perfect_p) {
      (*_perfect_p)[i] = b1;
      (*_perfect_p)[i+1] = b2;
    }

    // If we are generating errors, re-read the fragment using an
    // AlignmentProfile to create simulated reads with errors
    if ( _errgen_p ) {
      AlignmentProfile profile1, profile2;
      qualvector qual1, qual2;

      _errgen_p->GetErrorProfile(profile1, qual1, n);
      _errgen_p->GetErrorProfile(profile2, qual2, n);
      // Make sure the error profile isn't longer than the fragment!
      while ( profile1.GetRefBaseCount( ) > frag_length )
	_errgen_p->GetErrorProfile(profile1, qual1, n);
      while ( profile2.GetRefBaseCount( ) > frag_length )
	_errgen_p->GetErrorProfile(profile2, qual2, n);

      profile1.CreateRead  ( frag, 0, b1 );
      profile2.CreateRcRead( frag, frag.size( ) - 1, b2 );


      if ( jumping ) {
	b1.ReverseComplement( );
	b2.ReverseComplement( );
      }

      (*_qvv_p)[i]   = qual1;
      (*_qvv_p)[i+1] = qual2;
    }

    // we add to a local pair manager and append to the master one later in a mutex'd section

    (*_bvv_p)[i]   = b1;
    (*_bvv_p)[i+1] = b2;

  }

  if ( i_thread == 0 && _pm_p ) {
      for ( size_t i = 0; i < npairs; ++i ) {
	  _pm_p->addPair( 2*i, 2*i+1, N-2*n, static_cast<int>(round(dev)), "random", true);
      }
  }
}


// buildUnpairedReads -- allocates output structures and then calls
// CreateRandomReads() once per thread.

void ConstructionDReadGenerator::buildUnpairedReads( size_t read_size, size_t coverage )
{
  longlong genome_bases = 0;
  for ( size_t i = 0; i < _ref.size( ); i++ )
    genome_bases += _ref[i].size( );

  size_t nreads = int( round( float(genome_bases) * coverage / float(read_size) ) );

  _bvv_p->resize( _bvv_p->size() + nreads, BaseVec(read_size) );
  if ( _qvv_p ) _qvv_p->resize( _qvv_p->size() + nreads, QualVec(read_size) );
  if ( _perfect_p ) _perfect_p->resize( _perfect_p->size() + nreads, BaseVec(read_size) );
  if ( _refv_p ) _refv_p->resize( _refv_p->size() + nreads );


#pragma omp parallel for
  for ( size_t i_thread = 0; i_thread < _n_threads; ++i_thread ) {
    CreateRandomReads( nreads,  read_size, i_thread );
  }
}


// buildPairedReads and buildPairedJumps -- interface to buildPairs

void ConstructionDReadGenerator::buildPairedReads( size_t read_size, size_t coverage, size_t insert_size, double dev )
{
  this->buildPairs( read_size, coverage, insert_size, dev, -1, -1);
}


void ConstructionDReadGenerator::buildPairedJumps( size_t read_size, size_t coverage, size_t insert_size, double dev, int jmean, double jdev )
{
  this->buildPairs( read_size, coverage, insert_size, dev, jmean, jdev);
}


// buildPairs -- allocates output structures and then calls
// CreateRandomPairs() once per thread.  The latter will make reads if
// jmean = jdev = -1, otherwise jumps.

void ConstructionDReadGenerator::buildPairs( size_t read_size, size_t coverage, size_t insert_size, double dev, int jmean, double jdev )
{
  longlong genome_bases = 0;
  for ( size_t i = 0; i < _ref.size( ); i++ )
    genome_bases += _ref[i].size( );

  size_t npairs = size_t( round( float(genome_bases) * coverage / float( 2 * read_size ) ) );

  size_t nreads = 2*npairs;

  _bvv_p->resize( _bvv_p->size() + nreads, BaseVec(read_size) );
  if ( _qvv_p ) _qvv_p->resize( _qvv_p->size() + nreads, QualVec(read_size) );
  if ( _perfect_p ) _perfect_p->resize( _perfect_p->size() + nreads, BaseVec(read_size) );
  if ( _refv_p ) _refv_p->resize( _refv_p->size() + nreads );

  ForceAssert(_n_threads);	// must be non-zero


  cout << "Creating " << npairs << " read pairs for inserts of length "
       << insert_size << " +/- " << dev << " ";

  if ( jmean != -1 && jdev != -1 ){
    cout << "Jumping inserts are divided into fragments of length "
	 << jmean << " +/- " << jdev << "." << endl;
  }

  cout << endl;

#pragma omp parallel for
  for ( size_t i_thread = 0; i_thread < _n_threads; ++i_thread ) {
    CreateRandomPairs( npairs, insert_size, read_size, dev, jmean, jdev, i_thread );
  }

}



