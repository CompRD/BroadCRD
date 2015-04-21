/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// HowUnique.  For a given genome, compute the fraction of perfect N-base reads
// that can be aligned uniquely to it and the fraction of the genome that would
// be covered by them.  Files GENOME.fastb and GENOME.lookup must be provided.
//
// If a read has a false placement having <= D mismatches, the read is counted
// as not aligning uniquely.  Only alignments subsuming a K-base perfect match
// are seen.
//
// If SAMPLE is specified, only use that many reads.  In that case, coverage
// fraction can't be computed.
//
// If SAMPLE > 0, reads that involve ambiguous bases are skipped.  
// If SAMPLE == 0, ambiguous bases are not handled correctly.
//
// If COUNT_PERFECT=True, count mean number of perfect placements, and don't do
// anything else.  For this, SAMPLE must also be specified.  If also
// COUNT_PERFECT_VERBOSE=1, show number of placements of each read.  If 
// COUNT_PERFECT_VERBOSE=2, also show sequence.
//
// COUNT_PERFECT_FALSE: modifier for COUNT_PERFECT, to subtract one from each count
//
// We allow N to be a ParseIntSet-style list.

#include "Basevector.h"
#include "FeudalMimic.h"
#include "Intvector.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "lookup/ImperfectLookup.h"
#include "lookup/LookAlign.h"
#include "lookup/PerfectCount.h"
#include "lookup/PerfectLookup.h"
//#include "math/Functions.h"
#include "lookup/AlignCollector.h"
#include "lookup/HitReceiver.h"
#include "random/MersenneTwister.h"
#include "math/SparseArray.h"
#include "bias/BiasGC.h"


// Below copied from ShortQueryLookup.cc
template<class AlignCollectorType>
void LookupQueries(AlignCollectorType &aligns, lookup_table &look,
                   const vecbasevector &seqs,
                   unsigned int MAX_INDEL_LEN, unsigned int MAX_FREQ,
                   int nPasses, int AT_START)
{ 
  GlobalHitReceiver<AlignCollectorType> 
    receiver(look, seqs, &aligns, MAX_INDEL_LEN, AT_START);
  look.FindHits(seqs, receiver, MAX_FREQ, nPasses);
}


const char* DOC =
"For a given genome, generates N-base synthetic reads, aligns them to this genome, and "
"computes some alignment-related statistics: average number of perfect matches per read, "
"or the fraction of uniquely aligned reads and the fraction of the genome that would be "
"covered by them, or per-base coverage of the genome by these reads, or the best and "
"next-best alignment mismatch counts for these reads. The source sequence from which the "
"reads are generated can be either the genome itself (default), or a separate source "
"sequence fastb file (if specified). The reads can be either sampled randomly from the "
"source sequence, or come from every position on the source sequence. Finally, generated "
"synthetic reads can be randomly mutated to contain fixed number of mismatches";



/// Generates random collection of <sample_size> reads of length <N> each and
/// \em adds them to the <reads> collection. Random
/// positions on the specified reference <source_seq> are selected and 
/// subsequences of length <N> are read off those positions. If <ambig> is not empty,
/// all the returned subsequences are guaranteed not to contain any ambiguous bases.
/// If trim_left and/or trim_right are specified, the returned sequences are guaranteed
/// not to overlap with first trim_left and last trim_right bases on every contig
/// (sequence) from <source_seq>. Generated <sample_size> sequences are returned in <reads>.
/// If GC counts container is provided, GC counts of W-basepair windows anchored at every generated
/// read position will \em added to it; order of gc counts added to <gc_cnts> corresponds to the
/// order reads are generated and added. If the read is too close to the contig end, so that 
/// window of length W can not be placed, the corresponding gc count value will be set to -1.

void GenerateRandomReadSample( 
        unsigned int sample_size, ///< generate that many reads
        unsigned int N, /// generate reads of this length
        vecbasevector & reads, ///< place generated reads into here
        vecbasevector & readsB, ///< place paired generated reads into here
        vec<int> & insert_sizes,
	const vecbasevector & source_seqs, ///< generate reads from random positions on the source
	const vecbitvector & ambig, ///< ambiguous bases in source_seqs (can be empty)
	unsigned int trim_left ,  ///< no reads should start within first trim_left bases on each contig
	unsigned int trim_right, ///< no reads should end within last trim_right bases on each contig
	// unsigned int trim_left = 0 ,  ///< no reads should start within first trim_left bases on each contig
	// unsigned int trim_right = 0, ///< no reads should end within last trim_right bases on each contig
        Bool paired,
        double target_insert_size, double insert_stddev, int insert_max,
	vec<int> * gc_cnts = 0, ///< gc counts for each generated read (synchronized with reads)
	int W=50 ///< window size for gc content computation
	) {

	     vec<unsigned int> contig_start;
	     contig_start.push_back(0);
	     for ( size_t i = 0; i < source_seqs.size( ); i++ ) {
	       contig_start.push_back( contig_start.back( ) + source_seqs[i].size( ) - trim_left - trim_right );
	     }

	     for ( unsigned int u = 0; u < sample_size; u++ ) {    
  	         // generate random position into the whole (concatenated) genome:
	         unsigned int r = genrand64_int64() % contig_start.back( ); 
		 size_t g;
		 // find the contig the position r falls into:
		 for ( g = 0; g < source_seqs.size( ); g++ ) {
		     if ( r < contig_start[g+1] ) break;
		 }
		 if ( r + N > contig_start[g+1] ) {    
		     // if the end of the read placed at position r falls over the
		     // contig end, then we do not want it; decrement read counter and retry:
		     --u;
		     continue;    
		 }

                 int current_insert_size;

                 if ( paired ) {
                     // int OK, since insert size should never be that long.
                     current_insert_size = (int)roundl(genrand64_Box_Mueller_Gaussian(target_insert_size, insert_stddev));
                     if (current_insert_size > insert_max)
                        current_insert_size = insert_max;
                     if ( r + 2*N + current_insert_size > contig_start[g+1] || 
                          r + 2*N + current_insert_size < contig_start[g] ) {

                        // If 2nd half of pair falls off the end or goes
                        // before the beginning of the contig, try again.
                        --u;
                        continue;
                     }
                 }

		 if ( ! ambig.empty() ) {
		     // if ambiguous bases are specified for the reference genome:
                     bool hit_ambig = false;
		     // find out if our random read hits an ambiguous base:
		     for ( unsigned int b = 0; b < N; ++b ) {
                        // if in paired mode the other half has an 
                        // ambiguity in it, reject.
                        if ( ambig[g][r-contig_start[g]+trim_left+b] ||
                             (paired && ambig[g][r-contig_start[g]+trim_left+b+current_insert_size]) ) {
                          hit_ambig = true;
                          break;
                        }
		     }

		     if ( hit_ambig ) {
		        // yes, the read we just generated contains an ambiguous base, we 
		        // do not want it; decrement read counter and retry 
                        --u;
                        continue;
		     }
		 }

		 // now we are satisfied with random read we just generated; save it

		 unsigned int pos = r-contig_start[g]+trim_left; // position on contig
		 if ( gc_cnts != 0 ) { // if we have to save gc counts too:
		   if ( pos + W > source_seqs[g].size() ) { // too short to place window for gc counting
		     gc_cnts->push_back(-1); // gc_cnts must be synchronized with reads; so we
                                             // store "NA" value for current read
		   } else {
		     gc_cnts->push_back( source_seqs[g].GcBases(pos, pos+W) );
		   }
		 }

		 static basevector b; 
		 b.SetToSubOf( source_seqs[g], r - contig_start[g] + trim_left, N );
		 reads.push_back_reserve(b);    
                 if ( paired ) {
		    static basevector bB; 
		    bB.SetToSubOf( source_seqs[g], r - contig_start[g] + trim_left + current_insert_size, N );
		    readsB.push_back_reserve(bB);    
                    insert_sizes.push_back(current_insert_size);
                 }

	     } // end of for ( int u = 0 ; u < SAMPLE ; u++ ) - 
	       // loop that generates SAMPLE random reads    
}

/// Introduce exactly num_errors distinct random "single" (can accidentaly land next 
/// to each other!) mutations into the passed nucleotide sequence.
void MutateSequence(basevector & seq, unsigned int num_errors) {
  if ( num_errors == 0 ) return; // paranoidal sanity check in case meaningless value is passed
    vec<unsigned int> pos;
    pos.reserve(num_errors);
    unsigned int seq_size = seq.size();
    // generate exactly num_errors distinct random positions in the read:
    while(1) {    
        // generate one random position
        unsigned int p = genrand64_int64() % seq_size;
	// oops, we already have error at this position, retry:
	if ( Member( pos, p ) ) continue; 
	
	pos.push_back(p);
	if ( pos.size( ) == num_errors ) break;    // done!
    }
    // add a random number from [1,3] interval to each base at
    // the error positions we just generated:
    for ( unsigned int k = 0; k < num_errors; k++ ) {    
        unsigned int p = pos[k];
	unsigned int add = ( genrand64_int64( ) % 3 ) + 1; 
	seq.Set( p, ( seq[p] + add ) % 4 );    
    }    
}


template <class AlignCollector>
void ComputeCoverage(const AlignCollector & aligns, VecUIntVec & counts, bool unique_only = true) {
  for ( unsigned int i = 0; i < aligns.size( ); i++ ) {    
    if ( unique_only && ! aligns.UniquelyAligned(i) ) continue;
    const look_align& la = aligns.Align(i);
    int start = la.StartOnTarget( ), stop = la.EndOnTarget( );
    int tig = la.target_id;
    for ( int j = start; j < stop; j++ ) {
      ++counts[tig][j];
    }
  }
}

template <class AlignCollector>
void ComputeCoverage(const AlignCollector & aligns, vecbitvector & flags, bool unique_only = true) {
  for ( unsigned int i = 0; i < aligns.size( ); i++ ) {    
    if ( unique_only && ! aligns.UniquelyAligned(i) ) continue;
    const look_align& la = aligns.Align(i);
    int start = la.StartOnTarget( ), stop = la.EndOnTarget( );
    int tig = la.target_id;
    for ( int j = start; j < stop; j++ ) {
      flags[tig].Set(j,true);
    }
  }
}

/// translates absolute coordinate (offset) into contig:position coordinate on an arbitrary 2-dimensional
/// array structure <seqs>. Container <seqs> must conform to the following contract:
///   - seqs.size() is defined
///   - seqs.operator[] is defined (<contig> coordinate will be sought along this direction)
///   - seqs[i].size() is defined for each contig <i>
template <class Container2D>
void AbsPositionToContigPos(longlong abs_pos, unsigned int & contig, unsigned int & pos, const Container2D & seqs ) {
  longlong j = 0;
 
  for ( contig = 0 ; contig < (unsigned int) seqs.size() ; contig++ ) {
    if ( j + seqs[contig].size()  > abs_pos ) {
      pos = (unsigned int) (abs_pos - j);
      break;
    } 
    j += seqs[contig].size();
  }
}

///
template <class AlignCollector>
void ComputeAlignability(
 		     vecbitvector & flags,
		     const AlignCollector & a, ///<aligns
		     int N,  ///<read length
		     ulonglong start, /// start position on the read source reference
		     int CHUNK, ///< CHUNK size
		     ulonglong end_position ///< absolute end position on source_seqs 
		     ) {
  unsigned int contig = 0 ;
  unsigned int pos = 0;
  unsigned int nreads = 0;

  AbsPositionToContigPos(start, contig, pos, flags);

  // start generating reads from the current contig:pos and on
  for ( ; contig <  static_cast<unsigned int>(flags.size( )) && 
	  nreads < static_cast<unsigned int>(CHUNK) &&
	  start < end_position ;
	contig++, pos=0 ) {    

    // we get inside the next for loop only if we can generate a read of length N off position pos:

    for (  ; ( pos <= static_cast<unsigned int>(flags[contig].size( )) - static_cast<unsigned int>(N) ) &&
	              start < end_position &&
	     nreads < static_cast<unsigned int>(CHUNK)  ; 
              pos++, nreads++, start++ ) {    
      // hooray, we could generate a read from poisition pos earlier, let's see if it aligns!
      flags[contig].Set(pos, a.UniquelyAligned(nreads)); 
    }    
    if ( flags[contig].size() - pos < static_cast<unsigned int>(N) ) { 
      // not enough space from pos to the end of contig to generate a single read:
      // just mark all positions till the end of the contig as unalignable and advance to the next contig 
      for ( ; pos < flags[contig].size() ; pos++, start++ ) {
	flags[contig].Set(pos,False);
      }
      continue;
    }

		 
  }   
}



/// Generates a set of at most <CHUNK> reads of length <N> starting at every base in the reference <source_seqs> from
/// base with absolute offset <current_position> (inclusive) up to base with absolute offset <end_position>.
/// Will stop and return set of reads accumulated so far if <end_position> is encountered before <CHUNK> reads are generated.
/// Technically, it is more accurate to say that an <attempt> is made to generate a read from each and every position on 
/// the reference in the interval of absolute offsets [current_position, end_position). If it is impossible to generate a read
/// of length <N> from a given position (too close to the contig boundary), the position will be skipped and no read will be generated.
/// This method modifies <current_position>: upon return, <current_position> points to the first absolute position on the reference
/// from which next read should be generated. Note: this may be a position from which actual next read should and *can* be generated;
/// if, for instance, we got <CHUNK> reads and just reached (contig.size - N) position on the contig (no more reads can be generated
/// from this contig till its boundary), <current_position> upon return can be set to the beginning of the next contig right away.

void GenerateChunkOfReads(int N,  ///<read length
		     int CHUNK, ///< generate at most CHUNK reads (or less if end_position is encountered first
		     vecbasevector & reads, ///< generated reads will be <added> into this container
		     const vecbasevector & source_seqs, ///< set of sequences (contigs) to generate reads from
		     const vecbitvector & ambig, ///<ambiguous bases in source_seqs (can be empty)
		     longlong & current_position, ///< absolute starting position (i.e. first read to generate) on source_seqs
		     longlong end_position ///< absolute end position on source_seqs (no reads will be generated starting at or past this position)
		     ) {
  int nreads = 0; // counts reads generated so far during current invocation
  unsigned int contig = 0 ;
  unsigned int pos = 0;
  int n_ambig_reads = 0; // number of reads with ambiguous bases

  AbsPositionToContigPos(current_position, contig, pos, source_seqs);

  cout << "from " << contig << ":" << pos ;
  unsigned int last_c(0), last_p(0); // just for the fun of it, keep and print later last pos. used
  // start generating reads from the current contig:pos and on
  for ( ; contig < (unsigned int) source_seqs.size( ) && 
	  nreads < CHUNK &&
	  current_position < end_position ;
	contig++, pos=0 ) {    

    last_c = contig;

    // we get here only if we can generate a read of length N off position pos:

    for (  ; ( pos <= source_seqs[contig].size( ) - N) &&
	              current_position < end_position &&
	              nreads < CHUNK  ; 
              pos++, nreads++, current_position++ ) {    
      static basevector b;
      b.SetToSubOf( source_seqs[contig], pos, N ); // extract read
      reads.push_back_reserve(b);
      if ( ! ambig.empty() ) {
	  for ( int z = 0 ; z < N ; z++ ) {
	      if ( ambig[contig][pos+z] ) {
	          n_ambig_reads++;
	          break;
	      }
	  }
      }
    }    
    last_p = pos - 1; // save last position, from which a read was generated

    if ( source_seqs[contig].size() - pos < static_cast<unsigned int>(N) ) { 
      // not enough space from pos to the end of contig to generate a single read:
      // just advance current_position to the next contig
      current_position+=(source_seqs[contig].isize()-pos);
    }

		 
  }   

  cout << " to " << last_c << ":" << last_p << endl;
  cout << nreads << " reads generated (" << n_ambig_reads << " contain ambiguous bases)" <<  endl;
  flush(cout);
}


int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandDoc( DOC );
     CommandArgument_String_Doc(GENOME,"Reference lookup table (no .lookup extension) against which generated "
				"reads will be aligned");
     CommandArgument_Int_Doc(D, "alignments are deemed unique if next-best has at least D+1 more errors");
     CommandArgument_String_Doc(N, "length of reads");

     CommandArgument_String_OrDefault_Doc(GENOME2, "","Additional reference lookup table. This is a HACK. "
					   "Generated reads will be aligned against *both* GENOME *and* GENOME2 "
					   "and uniqueness of alignments or numbers of perfect matches will be "
					   "counted against *all* such alignments regardless of what reference the fell onto");
     CommandArgument_String_OrDefault_Doc(FASTB, "", 
	    "Collection of sequences, from which the reads will be generated. "
            "If not specified then genome itself, GENOME.fastb (all contigs), will be used." );
     CommandArgument_String_OrDefault_Doc(MODIFY,"","Apply specified modification (e.g. BS) to all generated reads, "
					  "prior to introducing errors (if any) and attempting to align") 
     CommandArgument_Int_OrDefault_Doc(K, 12, "matches must subsume a perfect K-mer");
     CommandArgument_Int_OrDefault_Doc(SAMPLE, 0, "how many samples? if 0, sample ALL N-mers");
     CommandArgument_String_OrDefault_Doc(QUEUE,"","If specified, alignment jobs will be broken into chunks and sent to the LSF queue");
     CommandArgument_Int_OrDefault_Doc(LSF_CHUNK,200000,"Submit this many reads per job when submitting alignment job to lsf queue in chunks");
     CommandArgument_Bool_OrDefault_Doc(COUNT_PERFECT, False,
                      "instead of computing the # of uniquely placed reads, show average # of perfect matches per read");
     CommandArgument_Bool_OrDefault_Doc(COUNT_PERFECT_FALSE, False,
                                        "show average # of perfect but false matches per read");
     CommandArgument_Int_OrDefault_Doc(COUNT_PERFECT_VERBOSE, 0,"print to stdout: 0 - nothing, "
				       "1 - perfect placement count for each read (one column) "
				       "2 - fasta with actual read sequences and with placement counts shown in the '>' lines");
     CommandArgument_Bool_OrDefault_Doc(WRITE,False,"Write all additional results into files (both _COV and _DIFF)");
     CommandArgument_Bool_OrDefault_Doc(WRITE_COV,False,"Write full coverage counts into a file; coverage is computed on GENOME");
     CommandArgument_Bool_OrDefault_Doc(WRITE_DIFF,False,"Write distances to next-best hit into a file");
     CommandArgument_Bool_OrDefault_Doc(WRITE_ALIGNABILITY,False,"Write read alignability bit flags; the bitflag array has "
					"the same structure as the source, from which the reads are generated: either GENOME or FASTB "
					"*after* trimming with TRIM_LEFT, TRIM_RIGHT, if any");
     CommandArgument_Bool_OrDefault_Doc(GC_STATS,False,"Compute GC count histogram");
     CommandArgument_Int_OrDefault_Doc(W,50,"Window size for GC stats computations");
     CommandArgument_String_OrDefault_Doc(OUT_HEAD,"","If specified, filename head for output file(s); "
					  "otherwise will use FASTB file head (if given) or GENOME");
     CommandArgument_Int_OrDefault_Doc(CHUNK,infinitely_many,"Generate reads by chunks of this size");
     CommandArgument_Int_OrDefault_Doc(ERRORS, 0,
          "Introduce ERRORS mismatches into every read.  In paired mode, introduce ERRORS mistmatches in each half of every read.");
     CommandArgument_Int_OrDefault_Doc(PRECISION, 3, "for printing results");
     CommandArgument_Int_OrDefault_Doc(SEED,0,"random number generator seed; use same seed "
				       "to reproduce \"random\" results");
     CommandArgument_Int_OrDefault_Doc(TRIM_LEFT,0,"Trim bases from the left side of each source "
				       "sequence/contig before generating reads");
     CommandArgument_Int_OrDefault_Doc(TRIM_RIGHT,0,"Trim bases from the right side of each source "
				       "sequence/contig before generating reads");
     CommandArgument_String_OrDefault_Doc(START,"0","Generate/sample reads from [START,END) interval on the source (0 based, continuous enumeration across all trimmed contigs)");
     CommandArgument_String_OrDefault_Doc(END,"-1","default is whole genome");
     CommandArgument_String_OrDefault_Doc(MAP,"","Map file with list of genomic regions <contig> <start> <stop>, "
                                                 "one per line, to generate reads from");
     CommandArgument_Bool_OrDefault_Doc(FWRC,True,"Align against GENOME (and GENOME2 if specified) in both directions. "
					  "If false, only forward alignments agains GENOME(+GENOME2) will be sought");
     CommandArgument_Bool_OrDefault(PAIRED, False);
     CommandArgument_Double_OrDefault_Doc(AVG_INSERT, 0.0, "Average insert size of paired reads");
     CommandArgument_Double_OrDefault_Doc(STDDEV_INSERT, 0.0, "Standard deviation of [randomly choosen] insert size of paired reads");
     CommandArgument_Int_OrDefault_Doc(MAX_INSERT, 1000000000, "Maximum length (all inserts longer than this becomes this size) of [randomly choosen] insert size of paired reads");
     CommandArgument_UnsignedInt_OrDefault_Doc(MAX_INDEL_LEN,0,
          "Don't look for alignments with "
          "indels larger than this.  Used to limit bandwidth "
          "for banded Smith-Waterman aligner.");
     CommandArgument_Int_OrDefault_Doc(MAX_FREQ,0,
          "if positive, only use kmers occurring at most "
          "this many times in the target genome.  HACK: if maxFreq is "
          "negative, interpret it as -B, where B is the blocksize.  "
          "For each consecutive blocks of kmers, starting at the "
          "beginning of the read, take the lowest (nonzero) frequency "
          "kmer in that chunk, and put only those kmers into the query "
          "set.  That is, take one kmer from each block "
          "[0,B), [B,2B), ... .  This is a much less sensible but "
          "sometimes much faster way to search for placements in "
          "highly repetitive genomes. ");
     CommandArgument_Int_OrDefault_Doc(INSERT_SLOP, 0, 
          "Second pair of reads may be off in alignment distance by this much and still be considered the same.");
     CommandArgument_Double_OrDefault_Doc(MIN_INSERT_STDDEV, 3.0, 
          "Second pair of reads must be at least this far away else not unique.");

     EndCommandArguments;

     if (PAIRED && SAMPLE == 0) {
             cout << "Paired mode does not run with SAMPLE = 0, due to long run-times.  Exiting.\n";
             exit(1);
     }
     if (PAIRED && GENOME2 != "") {
             cout << "Paired mode does not run with GENOME2 option.  Not implemented yet.  Exiting.\n";
             exit(1);
     }

     if (PAIRED && (AVG_INSERT == 0.0 || STDDEV_INSERT == 0.0 || MAX_INSERT == 0)) {
             cout << "Paired mode requires AVG_INSERT, STDDEV_INSERT, and MAX_INSERT to be set to something other than zero.  Exiting.\n";
             exit(1);
     }
     if (PAIRED && COUNT_PERFECT) {
             cout << "Counting perfect matches not implemented in paired mode.  Exiting.\n";
             exit(1);
     }

     if ( WRITE_ALIGNABILITY && SAMPLE != 0 ) {
       cout << "WRITE_ALIGNABILITY can not be used with partial sampling of the source sequence (SAMPLE=)" << SAMPLE << endl;
       exit(1);
     }

     vec<int> Ns;
     ParseIntSet( N, Ns );

     if (COUNT_PERFECT_FALSE) COUNT_PERFECT = True;
     if (COUNT_PERFECT_VERBOSE > 0 ) COUNT_PERFECT=True;
     //    if (COUNT_PERFECT) ForceAssertGt( SAMPLE, 0 );

     if ( COUNT_PERFECT && GC_STATS ) {
       cout << "GC statistics is currently computed only in unique "
	    << "alignment mode, turn off COUNT_PERFECT" << endl;
       exit(1);
     }

     if ( SAMPLE == 0 && GC_STATS ) {
       cout << "GC histograms are currently computed only in SAMPLE mode" << endl;
       exit(1);
     }

     if ( PAIRED && GC_STATS ) {
       cout << "GC_STATS currently not supported in the PAIRED mode" << endl;
       exit(1);
     }

     if ( ( WRITE || WRITE_COV || WRITE_DIFF || GC_STATS ) && OUT_HEAD == "" ) {
       if ( FASTB == "" ) {
	 OUT_HEAD = GENOME;
       } else {
	 OUT_HEAD = FASTB.SafeBefore(".fastb");
       }
     }

     if ( WRITE ) {
       WRITE_DIFF = True;
       WRITE_COV = True;
     }


     AlignDir align_direction = FW_OR_RC;
     if ( ! FWRC ) align_direction = FW;

     init_genrand64(SEED);

     unsigned int ref1_contigs = 0; // number of contigs in the GENOME reference

     // source to generate reads from:
     vecbasevector source_seqs( FASTB == "" ? GENOME + ".fastb" : FASTB );
     ref1_contigs = source_seqs.size();
     vecbasevector reference_seqs;
     if (FASTB != "") {
       reference_seqs.ReadAll(GENOME + ".fastb");
       ref1_contigs = reference_seqs.size(); // now *this* is the actual reference...
     }
     vecbasevector & ref = ( FASTB == "" ? source_seqs : reference_seqs );

     // ambiguous base information is used only when generating reads so it makes sense
     // to look for fastamb couterpart of the actual source sequences (FASTB or GENOME, whichever is requested)
     String ambig_file = ( FASTB=="" ? GENOME + ".fastamb" : FASTB.SafeBefore(".fastb")+".fastamb");
     vecbitvector ambig;
     if ( IsRegularFile( ambig_file ) ) {
       cout << "Ambiguous base data is available" << endl;
       ambig.ReadAll( ambig_file );
     }

     if ( TRIM_LEFT != 0 || TRIM_RIGHT != 0 ) {
       if ( PAIRED ) {
               cout << "Warning: code has not been tested with TRIM != 0.\n.";
       }
       basevector b;
       bitvector v;
       for ( unsigned int i = 0 ; i < (unsigned int)source_seqs.size() ; i++ ) {
	 if ( source_seqs[i].isize() < TRIM_LEFT +TRIM_RIGHT ) {
	   source_seqs[i].resize(0); 
	   if ( !ambig.empty() ) ambig[i].resize(0);
	 } else {
	   b.SetToSubOf(source_seqs[i],TRIM_LEFT, source_seqs[i].size()-TRIM_LEFT-TRIM_RIGHT);
	   source_seqs[i] = b;
	   if ( ! ambig.empty() ) {
	     v.SetToSubOf(ambig[i],TRIM_LEFT, ambig[i].size()-TRIM_LEFT-TRIM_RIGHT);
	     ambig[i] = v;
	   }
	 }
       }
       cout << "Source sequences trimmed" << endl;
       flush(cout);
     }

     longlong total_source_bases = 0;
     for ( unsigned int i = 0 ; i < ((unsigned int) source_seqs.size()) ; i++ ) {
         total_source_bases += source_seqs[i].size() ;
     }  
     if ( total_source_bases == 0 ) {
       cout << "No bases found in source sequences. If you are using TRIM_ make sure you are not trimming too much!" << endl;
       exit(0);
     }

     if ( ! START.IsInt() || ! END.IsInt() ) {
       cout << "START, END must be specified as integer numbers (K, M modifiers allowed)" << endl;
       exit(1);
     }

     longlong end_on_source= END.Int();
     if ( end_on_source == -1 || end_on_source > total_source_bases ) {
       end_on_source = total_source_bases;
     }
     
     longlong start_on_source = START.Int();
     if ( start_on_source >= end_on_source ) {
       cout << "START position exceeds END position: nothing to do" << endl;
       exit(0);
     }


     for ( int in = 0; in < Ns.isize( ); in++ )  {    
         int N = Ns[in];
         // N is length of read

	 size_t num_reads = 0;
	 if ( SAMPLE > 0 ) num_reads = SAMPLE;
	 else {
	   unsigned int contig1, pos1=0;
	   AbsPositionToContigPos(start_on_source, contig1, pos1, source_seqs);
	   unsigned int contig2, pos2=0;
	   AbsPositionToContigPos(end_on_source-1, contig2, pos2, source_seqs);

	   if ( contig1 != contig2 ) {
	     if (  source_seqs[contig1].size() - pos1 >= (unsigned int)N ) {
	       num_reads += source_seqs[contig1].size() - pos1 - N + 1;
	     }
	     if (  source_seqs[contig2].size() >= (unsigned int)N ) {
	       num_reads += Min(pos2, source_seqs[contig2].size()-N) + 1 ;
	     }
	     for ( contig1++ ; contig1 < contig2 ; contig1++ ) {
	       if ( source_seqs[contig1].size() < (unsigned int)N ) continue;
	       num_reads += source_seqs[contig1].size() - N + 1;
	     }  
	   } else {
	     if ( source_seqs[contig1].size() >= (unsigned int)N ) {
	       num_reads = Min(pos2, source_seqs[contig1].size()-N) -
		 Min(pos1, source_seqs[contig1].size()-N) + 1 ;
	     }
	   }
	   cout << "Total of " << num_reads << " reads will be generated" << endl;
	   if ( num_reads == 0 ) {
	     cout << "With the parameters specified no reads can be generated, aborting" << endl;
	     exit(0);
	   }
	 }
	 // num_reads now holds total number of reads that will be generated

	 // yes/no flags showing if any particular base is covered;
	 // we predefine it here so that results from multiple
	 // chunks could be accumulated
	 vecbitvector * coverage_flags=0; 
	 if ( ! COUNT_PERFECT ) {
	   // if COUNT_PERFECT is set, we are not going to 
	   // compute coverage!
	   coverage_flags = new vecbitvector();
	   if (FASTB == "")
	     Mimic( source_seqs, *coverage_flags );
	   else 
	     Mimic( ref, *coverage_flags );
	 }
	 
	 // actual coverage count; we define it here so that
	 // results from multiple chunks could be accumulated
	 // but actually allocate/resize only if it is going to be used
	 VecUIntVec * full_coverage=0;
	 if ( WRITE_COV ) {
	   // actual coverage counts for each base; used only if WRITE_COV=true
	   full_coverage = new VecUIntVec();
           if (FASTB == "")
	      Mimic ( source_seqs, *full_coverage );
           else
	      Mimic ( ref, *full_coverage );
	 }
	 vecbitvector alignability;
	 if ( WRITE_ALIGNABILITY ) {
	   Mimic( source_seqs, alignability ) ;
	 }


	 // open stream if writing alignability counts is requested:
	 ofstream *diffout = 0;
	 if ( WRITE_DIFF ) {
	   diffout = new ofstream();
	   OpenOfstream( *diffout, "diffout", OUT_HEAD + "." +ToString(N)+".errDiff" );
	 }

	 // next variable is used to remember the position on the
	 // reference, from which reads are generated: if the reference is too 
	 // big, then we may want to generate reads in chunks. In that case
	 // we should be able to generate next chunk of reads starting from
	 // the stored position. we store the position as one number running
	 // continously over all contigs
	 longlong current_position = start_on_source;

	 // used to keep the continuous read niumber across all chunks; used only for
	 // WRITE_DIFF and WRITE_ALIGNABILITY
	 unsigned int read_number = 0;


	 double perfect_placements = 0.0;
	 size_t multiple_perfect = 0;
	 longlong global_multiple_max = 0;
	 size_t total_aligned = 0;

	 // these will be used to count funky placements of reads with errors
	 // (unused if ERRORS==0):
	 ulonglong total_better = 0; 
	 ulonglong total_worse = 0;
	 ulonglong unique_better = 0;
	 ulonglong unique_worse = 0;

	 unsigned int chunk_no = 0;


	 vec<int> GC_hist(W+1,0);
	 vec<int> GC_hist_aligned(W+1,0);


	 // until we generated reads from each and every position on the ref...
	 while ( current_position < end_on_source ) {
	   
	   // ##########################################################
	   // Generate chunk of reads from the sequence collection or SAMPLE random reads:

	   vecbasevector reads, readsB;
           vec<int> insert_sizes;
	   chunk_no++;

	   ulonglong chunk_start = current_position;

	   vec<int> * p_gc_counts = 0;
	   if ( GC_STATS ) {
	     p_gc_counts = new vec<int>();
	   }

	   if ( SAMPLE == 0 )  {    
	       // no sampling; generate next chunk of reads:
	     GenerateChunkOfReads(N, CHUNK, reads, source_seqs, ambig, current_position, end_on_source );
	     PRINT2(current_position, end_on_source);
	       
	   } else { 

	       // sampling requested; generate a random sample of reads with specified sample size 
	       GenerateRandomReadSample(SAMPLE,N,
					reads,readsB,
					insert_sizes,source_seqs,ambig,
					0 /*TRIM_LEFT*/,0 /*TRIM_RIGHT*/, // source_seqs are already trimmed!
					PAIRED,AVG_INSERT,
					STDDEV_INSERT,MAX_INSERT,p_gc_counts,W);
	       current_position = end_on_source ; // signal that no more chunks are needed
	       cout << SAMPLE << " reads randomly sampled from the genome." << endl;
	       if ( GC_STATS ) {
		 ForceAssertEq(p_gc_counts->size(), (unsigned int)reads.size());
		 for ( unsigned int i = 0 ; i < p_gc_counts->size() ; i++ ) {
		   if ( (*p_gc_counts)[i] == -1 ) continue; // gc count could not be computed
		   ++GC_hist[ (*p_gc_counts)[i] ]; // count the occurence of i-th gc count value
		 }
	       }
	   }


	   if ( MODIFY != "" ) {
	     if ( MODIFY == "BS" ) {
	       // emulate bisulfite treatment of the reads:
	       for ( unsigned int i = 0 ; i < (unsigned int)reads.size() ; i++ ) {
		 for ( unsigned int j = 0 ; j < reads[i].size() ; j++ ) {
		   if ( reads[i][j] == BASE_C ) reads[i].Set(j,BASE_T); // convert C to T
		 } 
		 cout << "All generated reads are bisulfite converted." << endl;
	       }
	       // that's all folks!
	     } else {
	       cout << "Unknown sequence modification " << MODIFY << " is requested. Exiting." << endl;
	       exit(1);
	     }
	   }
	   // introduce errors into generated reads
	   if ( ERRORS > 0 ) {
	     for ( unsigned int i = 0 ; i < (unsigned int)reads.size() ; i++ ) {
	       MutateSequence(reads[i],ERRORS);
               if ( PAIRED ) MutateSequence(readsB[i],ERRORS);
	     }
	     cout << ERRORS << " random substituions per read introduced" << endl;
	   }
	   // Reads generated
	   // ################################################
	   
	   // Now align the reads back and see what happens:
	   
	   
	   if (COUNT_PERFECT)  {    
	       // deal with counting perfect matches if asked:

	       vec<int> places;
               PerfectCountPlaces( reads, GENOME + ".lookup", align_direction , places );

	       int placemax = 0;
	       unsigned int m_perf = 0; // counter of reads with mult. placements from current chunk
	       for ( unsigned int i = 0 ;  i < places.size() ; i++ ) {
		 if ( places[i] > 1 ) m_perf++;
		 if ( places[i] > placemax ) placemax = places[i];
	       }

	       perfect_placements += Sum(places);

	       if ( reads.size() < num_reads ) {
		 cout << "<chunk_" << chunk_no << ">: ";
		 if ( GENOME2 != "" ) cout << "pass1: " ;
		 cout << ((double)Sum(places))/reads.size() 
		      << " perf. placements per read" << endl;
		 cout << "<chunk_" << chunk_no << ">: ";
		 if ( GENOME2 != "" ) cout << "pass1: " ;
		 cout << m_perf << " reads have multiple perf., placements (max. " 
		      << placemax << " placements)" << endl;
	       }  
	       if ( GENOME2 != "" ) {
		 vec<int> places2;
		 PerfectCountPlaces( reads, GENOME2 + ".lookup", align_direction , places2 );

		 int placemax2 = 0;
		 unsigned int m_perf2 = 0; // counter of reads with mult. placements from current chunk
		 for ( unsigned int i = 0 ;  i < places2.size() ; i++ ) {
		   if ( (places[i] + places2[i]) > 1 ) m_perf2++;
		   if ( (places[i]+places2[i]) > placemax2 ) placemax2 = (places[i]+places2[i]);
		 }

		 if ( reads.size() < num_reads ) {
		   cout << "<chunk_" << chunk_no << ">: pass2: ";
		   cout << ((double)Sum(places2)+Sum(places))/reads.size() 
		      << " perf. placements per read" << endl;
		   cout << "<chunk_" << chunk_no << ">: pass2: ";
		   cout << m_perf2 << " reads have multiple perf., placements (max. " 
		      << placemax2 << " placements)" << endl;
		 }  
		 multiple_perfect+=m_perf2; // total counter
		 if ( placemax2 > global_multiple_max) global_multiple_max = placemax2;

		 perfect_placements += Sum(places2);

	       } else {
		 multiple_perfect+=m_perf; // total counter
		 if ( placemax > global_multiple_max) global_multiple_max = placemax;
	       }
	       
	       if ( COUNT_PERFECT_VERBOSE == 1 ) {    
	         for ( unsigned int i = 0 ; i < places.size() ; i++ ) cout << places[i] << "\n";    
	       }
	       if ( COUNT_PERFECT_VERBOSE == 2 ) {    
		 for ( unsigned int i = 0 ; i < places.size() ; i++ )  {    
		   reads[i].Print( cout, "read_" + ToString(i) + " - occurs " 
				   + ToString(places[i]) + " times" );    
		 }    
	       }
	   } else {

	       // otherwise (if not COUNT_PERFECT) count all the ungapped 
	       // alignments (including imperfect ones):
	       //UniqueByErrDiffAlignCollector aligns(D,infinitely_many);

	       // We have to copy-paste same big chunk of code here only because
	       // we want different collectors here (one is nore efficient than another),
	       // and they are not polymorphic! This "design" should go!!!
	       if ( WRITE_DIFF ) { 
	           BestNextBestAlignCollector aligns(D,infinitely_many);

		   ImperfectLookup( GENOME + ".lookup", reads, aligns,align_direction);
		   unsigned int cnt = UniquelyAlignedCount(aligns);
		   if ( reads.size() < num_reads ) {
		     // printout partial stats only if we are doing chunks
		     cout << "<chunk_" << chunk_no << ">: ";
		     if ( GENOME2 != "" ) cout << "pass1: " ;
		     cout << cnt << " reads (" 
			  << PERCENT_RATIO( PRECISION, cnt,(unsigned int) reads.size( ) )
			  << " aligned uniquely" << endl;
		   }

		   if ( GENOME2 != "" ) {
		     BestNextBestAlignCollector aligns2(D,infinitely_many);
		     ImperfectLookup( GENOME2 + ".lookup", reads, aligns2,align_direction);
		     aligns2.AdjustTargetIds(ref1_contigs);
		     Combine(aligns,aligns2);

		     cnt = UniquelyAlignedCount(aligns);
		     if ( reads.size() < num_reads ) {
		       // printout partial stats only if we are doing chunks
		       cout << "<chunk_" << chunk_no << ">: pass2:";
		       cout << cnt << " reads (" 
			    << PERCENT_RATIO( PRECISION, cnt,(unsigned int) reads.size( ) )
			    << " aligned uniquely" << endl;
		     }
		   }

		   total_aligned += UniquelyAlignedCount(aligns);
		   
		   if ( ERRORS > 0 ) {
		     // if we place reads with errors, they may end at some
		     // funky locations: 1) a read generated from a given place with
		     // K mismatches may find itself a place elsewhere with fewer than
		     // K mismatches; 2) we may fail to find correct placement and get 
		     // a worse one - we want to count all these events here
		     for ( unsigned int i = 0 ; i < aligns.size() ; i++ ) {
		       if ( aligns.MinErrors(i) == -1  ) continue; // if we've seen an alignment at all
		       if ( aligns.MinErrors(i) < ERRORS ) {
			 // oops, we found a better alignment at wrong place! - 
			 // alignment at correct place should have had exactly ERROR mismatches!
			 total_better++;
			 if ( aligns.UniquelyAligned(i) ) unique_better++; // and it's unique!
		       } else {
			 if ( aligns.MinErrors(i) > ERRORS ) {
			   // oops, we failed to find correct alignment, and ended
			   // up with a worse one
			   total_worse++;
			   if ( aligns.UniquelyAligned(i) ) unique_worse++; // and it's unique!
			 }
		       }
		     }
		   }

		   // now accumulate alignments for the current chunk into coverage map
		   ComputeCoverage(aligns,*coverage_flags);
		   
		   if ( WRITE_COV ) {
		     // if we write actual coverage counts, it's time to accumulate
		     // these counts for alignments from the current chunk
		     ComputeCoverage(aligns,*full_coverage);
		   }
		   
		   if ( WRITE_DIFF ) {
		     // alignment best-next best mismatch difference is per read, not
		     // per base on the reference, so that we do not have to store
		     // this information and can dump it to disk immediately for every read
		     // from the current chunk:
		     for ( unsigned int i = 0 ; i < aligns.size() ; i++ ) {
		       (*diffout) << read_number << "\t" << aligns.MinErrors(i) 
				  << "\t" << aligns.NextErrors(i) << endl;
		       read_number++; // extrenal counter, makes sure we enumerate reads continuously across chunks
		     }
		   }
	       } else {
		   unsigned int cnt;
		   int MAX_ERR = 4;
	           UniqueByErrDiffAlignCollector aligns(D,MAX_ERR);

		   if ( QUEUE == "" ) { // don't use lsf, align directly:
		     ImperfectLookup( GENOME + ".lookup", reads, aligns,align_direction);
		   } else {
		     // spawn smaller lsf jobs, wait until completion:
		     String temp_reads_name = "tmp.howunique."+ToString(N)+".fastb";
		     reads.WriteAll(temp_reads_name);
		     AlignInChunksOnFarm_ILT(temp_reads_name,
					     GENOME + ".lookup",
					     "tmp.howunique",
					     "."+ToString(N)+".qltout",
					     QUEUE,
					     "ERR_DIFF=" + ToString(D) + " MAX_ERRS=" + ToString(MAX_ERR)+" WRITE_MINALIGNS=False" +
					     (MAX_FREQ > 0 ? (" MAX_FREQ="+ToString(MAX_FREQ)): ""), LSF_CHUNK
					     );
		     cout << "Alignment completed on LSF farm. Loading aligns..." ; flush(cout);
		     aligns.resize(reads.size());
		     LoadLookAligns("tmp.howunique."+ToString(N)+".qltout",aligns);
		     cout << "done" << endl; flush(cout);
		     //Remove(temp_reads_name);
		     //		     Remove("tmp.howunique."+ToString(N)+".qltout");
		   }

                   // Paired data structure below...
                   vecbasevector uniqueAcheckB, checkB, uniqueBcheckA;
                   vec<int> uniqueAcheckB_As_id, checkBorgId, uniqueBcheckA_Bs_id, uniqueBcheckA_Borg_id;
	           UniqueByErrDiffAlignCollector alignsB(D,ERRORS+3);
	           UniqueByErrDiffAlignCollector checkBalign(D,ERRORS+3);
	           UniqueByErrDiffAlignCollector alignsBA(D,ERRORS+3);

                   if ( !PAIRED ) {
                        cnt  = UniquelyAlignedCount(aligns);
                   } else {
                        alignsB.resize(reads.size());
                        // checkBalign will be resized by ImperfectLookup
                        alignsBA.resize(reads.size());
                        for (unsigned int id = 0; id < aligns.size(); ++id) {
                            if ( aligns.UniquelyAligned(id) ) {
                                // We'll check if there are corresponding
                                // B non-unique alignments
                                uniqueAcheckB.push_back_reserve(readsB[id]);
                                uniqueAcheckB_As_id.push_back(id);
                            } else {
                                // Does not align uniquely so delete alignment.
                                // We'll align the other side and see if they 
                                // align uniquely and if so, then we check A's
                                // non-unique alignments
                                aligns.Erase(id);
                                checkB.push_back_reserve(readsB[id]);
                                checkBorgId.push_back(id);
                            }
                        }
                        lookup_table look( GENOME + ".lookup" );
                        const int nPasses = 1;
                        const int AT_START = -1;
                        // First check if uniqueA's other pair read align 
                        // non-uniquely to the "right" places.
                        MaxErrDiffAlignCollector aligns_uniqueAcheckB(D,ERRORS+3);
                        LookupQueries(aligns_uniqueAcheckB, look, 
                                      uniqueAcheckB, MAX_INDEL_LEN, MAX_FREQ, 
                                      nPasses, AT_START);
                        aligns_uniqueAcheckB.Consolidate();
                        for (unsigned int i = 0; i < uniqueAcheckB_As_id.size(); i++) {
                            const look_align& la = aligns.Align(uniqueAcheckB_As_id[i]);
                            ho_interval exA = la.Extent2(), exB; 
                            const vec<look_align>& lb = aligns_uniqueAcheckB.Aligns(i);
                            bool good_pair = true;
                            unsigned int j, goodj = 0;
                            for (j = 0; j < lb.size(); ++j) {
                                exB = lb[j].Extent2();
                                int d = abs(Distance(exA, exB)) + N;
                                int d_diff = abs(d - insert_sizes[uniqueAcheckB_As_id[i]]);
                                if ( d_diff <= INSERT_SLOP ) 
                                    goodj = j;
                                else if ((d_diff > INSERT_SLOP) && 
                                         (d_diff <= MIN_INSERT_STDDEV*STDDEV_INSERT))
                                    good_pair = false;
                            }
                            if ( good_pair ) 
                                // j below is the last j that worked above.
                                alignsB.Insert(lb[goodj]);
                            else 
                                // Not a good pair as second half did not
                                // align in the right place, so delete 
                                // alignment in A
                                aligns.Erase(uniqueAcheckB_As_id[i]);
                        }

                        // Now, try to uniquely align all the alignments in 
                        // checkB
		        ImperfectLookup( GENOME + ".lookup", checkB, checkBalign,align_direction);
                        for (unsigned int id2 = 0; id2 < checkBalign.size(); ++id2) {
                            if ( checkBalign.UniquelyAligned(id2) ) {
                                uniqueBcheckA.push_back_reserve(reads[checkBorgId[id2]]);
                                uniqueBcheckA_Bs_id.push_back(id2);
                                uniqueBcheckA_Borg_id.push_back(checkBorgId[id2]);
                            } else {
                                // Neither end aligns uniquely :-(
                                checkBalign.Erase(id2);
                            }
                        }
                        MaxErrDiffAlignCollector aligns_uniqueBcheckA(D,ERRORS+3);
                        LookupQueries(aligns_uniqueBcheckA, look, 
                                      uniqueBcheckA, MAX_INDEL_LEN, MAX_FREQ, 
                                      nPasses, AT_START);
                        aligns_uniqueBcheckA.Consolidate();
                        for (unsigned int i2 = 0; i2 < uniqueBcheckA_Bs_id.size(); i2++) {
                            const look_align& lb = checkBalign.Align(uniqueBcheckA_Bs_id[i2]);
                            ho_interval exA, exB = lb.Extent2();
                            const vec<look_align>& la = aligns_uniqueBcheckA.Aligns(i2);
                            bool good_pair = true;
                            unsigned int j2, goodj2 = 0;
                            for (j2 = 0; j2 < la.size(); ++j2) {
                                exA = la[j2].Extent2();
                                int d = abs(Distance(exA, exB)) + N;
                                int d_diff = abs(d - insert_sizes[uniqueBcheckA_Borg_id[i2]]);
                                if (d_diff <= INSERT_SLOP) 
                                    goodj2 = j2;
                                else if ((d_diff > INSERT_SLOP) && 
                                         (d_diff <= MIN_INSERT_STDDEV*STDDEV_INSERT))
                                    good_pair = false;
                            }
                            if ( good_pair ) 
                                alignsBA.Insert(la[goodj2]);
                            else 
                                checkBalign.Erase(i2);
                        }
                        cnt = UniquelyAlignedCount(aligns);
                        cnt += UniquelyAlignedCount(alignsB);
                        cnt += UniquelyAlignedCount(checkBalign);
                        cnt += UniquelyAlignedCount(alignsBA);

                        // At the end, we have aligns, alignsB, 
                        // checkBalign, and alignsBA.
		   
                   } // END IF( ! PAIRED ) {} ELSE {

		   if ( reads.size() < num_reads ) {
		     // printout partial stats only if we are doing chunks
		     cout << "<chunk_" << chunk_no << ">: ";
		     if ( GENOME2 != "" ) cout << "pass1: " ;
		     cout << cnt << " reads (" 
			  << PERCENT_RATIO( PRECISION, cnt,(unsigned int) reads.size( ) )
			  << " aligned uniquely" << endl;
		   }

		   if ( GENOME2 != "" ) {
		     UniqueByErrDiffAlignCollector aligns2(D,ERRORS+3);
		     ImperfectLookup( GENOME2 + ".lookup", reads, aligns2,align_direction);

		     aligns2.AdjustTargetIds(ref1_contigs);
		     Combine(aligns,aligns2);

		     cnt=UniquelyAlignedCount(aligns);
		     if ( reads.size() < num_reads ) {
		       // printout partial stats only if we are doing chunks
		       cout << "<chunk_" << chunk_no << ">: pass2: ";
		       cout << cnt << " reads (" 
			    << PERCENT_RATIO( PRECISION, cnt,(unsigned int) reads.size( ) )
			  << " aligned uniquely" << endl;
		     }

		   } // if ( GENOME2 != "" )

		   if ( GC_STATS ) {
		     ForceAssertEq(aligns.size(), p_gc_counts->size());
		     for ( unsigned int i = 0 ; i < aligns.size() ; i++ ) {
		       if ( (*p_gc_counts)[i] == -1 ) continue; // se did not have count for i-th read
		       if ( ! aligns.UniquelyAligned(i) ) continue; // read is not uniquely aligned
		       ++GC_hist_aligned[ (*p_gc_counts)[i] ]; // read is uniquely aligned - count its 
		                                               // gc content into the appropriate bin
		     }
		   }

		   total_aligned += UniquelyAlignedCount(aligns);
                   if ( PAIRED ) {
                        total_aligned += UniquelyAlignedCount(alignsB);
                        total_aligned += UniquelyAlignedCount(checkBalign);
                        total_aligned += UniquelyAlignedCount(alignsBA);
                   }
		   
		   if ( ERRORS > 0 ) {
		     // if we place reads with errors, they may end at some
		     // funky locations: 1) a read generated from a given place with
		     // K mismatches may find itself a place elsewhere with fewer than
		     // K mismatches; 2) we may fail to find correct placement and get 
		     // a worse one - we want to count all these events here
		     for ( unsigned int i = 0 ; i < aligns.size() ; i++ ) {
		       if ( aligns.MinErrors(i) == -1  ) continue; // if we could not find any alignments at all
		       if ( aligns.MinErrors(i) < ERRORS ) {
			 // oops, we found a better alignment at wrong place! - 
			 // alignment at correct place should have had exactly ERROR mismatches!
			 total_better++;
			 if ( aligns.UniquelyAligned(i) ) unique_better++; // and it's unique!
		       } else {
			 if ( aligns.MinErrors(i) > ERRORS ) {
			   // oops, we failed to find correct alignment, and ended
			   // up with a worse one
			   total_worse++;
			   if ( aligns.UniquelyAligned(i) ) unique_worse++; // and it's unique!
			 }
		       }
		     }
                     if ( PAIRED ) {
		         for ( unsigned int i = 0 ; i < alignsB.size() ; i++ ) {
		           if ( alignsB.MinErrors(i) == -1  ) continue;
		           if ( alignsB.MinErrors(i) < ERRORS ) {
			     total_better++;
			     if ( alignsB.UniquelyAligned(i) ) unique_better++;
		           } else {
			     if ( alignsB.MinErrors(i) > ERRORS ) {
			       total_worse++;
			       if ( alignsB.UniquelyAligned(i) ) unique_worse++; 
			     }
		           }
		         }
		         for ( unsigned int i = 0 ; i < checkBalign.size() ; i++ ) {
		           if ( checkBalign.MinErrors(i) == -1  ) continue;
		           if ( checkBalign.MinErrors(i) < ERRORS ) {
			     total_better++;
			     if ( checkBalign.UniquelyAligned(i) ) unique_better++;
		           } else {
			     if ( checkBalign.MinErrors(i) > ERRORS ) {
			       total_worse++;
			       if ( checkBalign.UniquelyAligned(i) ) unique_worse++; 
			     }
		           }
		         }
		         for ( unsigned int i = 0 ; i < alignsBA.size() ; i++ ) {
		           if ( alignsBA.MinErrors(i) == -1  ) continue;
		           if ( alignsBA.MinErrors(i) < ERRORS ) {
			     total_better++;
			     if ( alignsBA.UniquelyAligned(i) ) unique_better++;
		           } else {
			     if ( alignsBA.MinErrors(i) > ERRORS ) {
			       total_worse++;
			       if ( alignsBA.UniquelyAligned(i) ) unique_worse++; 
			     }
		           }
		         }
                     }  // End of if ( PAIRED ) 
		   } // END of IF ( ERRORS > 0 )

		   // now accumulate alignments for the current chunk into coverage map
		   ComputeCoverage(aligns,*coverage_flags);
		   if ( WRITE_ALIGNABILITY ) {
		     ComputeAlignability(alignability, aligns, N, chunk_start, CHUNK, end_on_source);
		   }

                   if ( PAIRED ) {
		      ComputeCoverage(alignsB,*coverage_flags);
		      ComputeCoverage(checkBalign,*coverage_flags);
		      ComputeCoverage(alignsBA,*coverage_flags);
                   }
		   
		   if ( WRITE_COV ) {
		     // if we write actual coverage counts, it's time to accumulate
		     // these counts for alignments from the current chunk
		     ComputeCoverage(aligns,*full_coverage);
                     if ( PAIRED ) {
		         ComputeCoverage(alignsB,*full_coverage);
		         ComputeCoverage(checkBalign,*full_coverage);
		         ComputeCoverage(alignsBA,*full_coverage);
                     }
		   }
	       }
	   }
	   

	   if ( GC_STATS ) {
	     delete p_gc_counts;
	   }
	   // done with processing data for the current chunk, go get next chunk
	 } // end WHILE (! ALL_CHUNKS_DONE )
	 
	 // done with all chunks, now process the results:

	 cout << "all chunks processed" << endl;
	 flush(cout);

	 if ( COUNT_PERFECT ) {
	     // get the average number of perfect placements across all reads:
	     double meanp = perfect_placements/num_reads;
	     if (COUNT_PERFECT_FALSE) --meanp;

	     cout << "N = " << N << " --> ";
	     RightPrecisionOut( cout, meanp, PRECISION );
	     if ( COUNT_PERFECT_FALSE ) cout << " false ";
	     cout << " perfect placements per read" << endl;
	     cout << multiple_perfect << "(" 
		  << PERCENT_RATIO(PRECISION, multiple_perfect, num_reads) 
		  <<") reads have multiple perfect placements" << endl;
	     cout << "Max. " << global_multiple_max << " placements per read" << endl;
	 } else {

	     // if not counting perfect, then we are looking at coverage: 
             if ( !PAIRED ) {
	        cout << total_aligned << " of " << num_reads 
		     << " reads (" ;
	        cout << PERCENT_RATIO( PRECISION, total_aligned, num_reads )
		     << ") can be aligned uniquely" << endl;
             } else {
                // Paired has twice as many reads.  :-)
	        cout << total_aligned << " of " << 2*num_reads 
		     << " single-end reads (" ;
	        cout << PERCENT_RATIO( PRECISION, total_aligned, 2*num_reads )
		     << ") can be aligned uniquely" << endl;
             }
	     cout << setprecision(PRECISION) << 100.0 * Coverage(*coverage_flags)
		  << "% of genome covered by uniquely aligning reads" 
		  << endl;    
	     delete coverage_flags;

	     if ( WRITE_COV ) {
	         // if writing actual coverage counts to disk is requested:
	         Ofstream( covout, OUT_HEAD + "." + ToString(N) + ".cover.txt" );
		 for ( size_t t = 0 ; t < full_coverage->size() ; t++ ) {
		   for ( unsigned int p = 0 ; p < (*full_coverage)[t].size() ; p++ ) {
		     covout << t << "\t" << p << "\t" << (*full_coverage)[t][p] << endl;
		     }
		 }
		 covout.close();
		 delete full_coverage;
	     }	     

	     if ( WRITE_ALIGNABILITY ) {
	       alignability.WriteAll( OUT_HEAD +"." + ToString(N) + ".alignability.bits" );
	     }
	 
	     if ( WRITE_DIFF ) {
	       diffout->close();
	       delete diffout;
	     }
	 
	     if ( GC_STATS ) {
	       Ofstream( gc_out, OUT_HEAD +"." + ToString(N) +".gchist.txt" );
	       gc_out << "GC_content\treads\taligned_reads" << endl;
	       for ( int i = 0 ; i <= W ; i++ ) {
		 gc_out << i << "\t" << GC_hist[i] << "\t" << GC_hist_aligned[i] << endl;
	       }
	       gc_out.close();
	     }
	       
	     if ( ERRORS > 0 ) {
	       // if we generated reads with errors, we compute additional
	       // statistics for misplaced reads: a) reads that found better 
	       // places than those they were actually generated from; b) reads,
	       // for which we failed to find any placement as good as original one and
	       // could foind only a worse one (this may be due only to aligner heuristics
	       // of course, so it's somewhat less interesting).

	       cout << total_better << " reads (" 
		    << setprecision(PRECISION) << ((double)total_better*100.0)/num_reads << " %) have placements with < " 
		    << ERRORS << " mismatches" << endl;
	       cout << total_worse << " reads (" 
		    << setprecision(PRECISION) << ((double)total_worse*100.0)/num_reads << " %) have placements with > " 
		    << ERRORS << " mismatches" << endl;
	       cout << setprecision(PRECISION) 
		    << "Of those, "
		    <<((double)unique_better*100)/total_better << " % / "
		    <<((double)unique_worse*100)/total_worse 
		    << " % are also unique, respectively" << endl;
	     } // end if (ERRORS)

	 //	 } // end of if (SAMPLE==0)    
	 } // end of ELSE clause for if ( COUNT_PERFECT )
	 //	 delete coverage_flags;
	 //	 delete full_coverage;
     } // end of for ( int Ns = 0;...) loop over all requested read sizes    
}
