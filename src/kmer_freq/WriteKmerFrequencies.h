/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef __INCLUDE_kmer_freq_WriteKmerFrequencies_h
#define __INCLUDE_kmer_freq_WriteKmerFrequencies_h

#include "Basevector.h"
#include "CoreTools.h"
#include "kmers/KmerRecord.h"
#include "Bitvector.h"

/**
   File: WriteKmerFrequencies.h

   Write out the frequencies of all (or selected) kmers in the reads; see <WriteKmerFrequencies()>.

   @file
*/

/**
   Type Concept: KSelector

   A filter on the set of occurrences of a given kmer in the reads.
   Looks at the set of occurrences of a kmer in the reads, and tells <WriteKmerFrequencies()>
   whether to write out that kmer and its frequency or to ignore it.
   To be a model of KSelector, a class must have the members defined below.

   Type: record_t

   The type of kmer records to extract; the type must be a model of <SortKmersOutputRecord>.
   This type will depend on what information <operator()()> needs to decide whether to allow
   the given kmer.  For example, if we don't care about the kmer's occurrences in the reads,
   then the <kmer> class is enough, but if we do, we'll need <kmer_record>.

   Operator: operator()

   Tells whether the given kmer should be included in the output, given the set of
   occurrences of this kmer in the reads. Returns the frequency to record for the
   given kmer.
   
   > unsigned short operator() (const vec< record_t >& kmerOccurrences,
   >                            int fromKmerIdx, int toKmerIdx) const;

   Parameters:

      kmerOccurrences - the vector of records indicating kmer occurrences;
          the type of records must model the <SortKmersOutputRecord> type
	  concept.  only the records in the range [fromKmerIdx, toKmerIdx),
	  all representing different occurrences of the same kmer in the reads,
	  are used.

      fromKmerIdx, toKmerIdx - the range of records in kmerOccurrences that
          represent the set of occurrences of the given kmer in the reads.
	  the predicate operator should only look at those records.

   Returns:
      the frequency to record, or 0 if the kmer should not be included in the
      output of <WriteKmerFrequencies()>.
*/

// End: Section

/**
   Function: WriteKmerFrequencies
  
   Find all the k-mers in a given <.fastb> file and
   its reverse complement.  Output these kmers and their
   multiplicities, capped at 65535, optionally excluding kmers that
   appear only once or that are rejected by a <KSelector> predicate.

   Template parameters:

      KSHAPE - the <shape of kmers> to extract from the reads
      KSELECTOR - the type of the kselector argument, which will
        select the kmers whose frequencies will be included in the output;
	must be a model of <KSelector>.
  
   Parameters:
  
      seqs - the reads or genome parts in which we want to gather kmers and their
         frequencies
      filename - a <KmerShortMap file> containing (kmer, frequency) pairs for each kmer
      include_unique - whether kmers appearing only once in _seqs_ should be included
        in the output.
      kselector - a filter on kmers, that looks at the set of occurrences
        of a given kmer and decides whether to include that kmer in the output.
*/
template <class KSHAPE, class KSELECTOR>
void
WriteKmerFrequencies( const vecbasevector& seqs,  const String& filename,
                      const bool include_unique, const KSELECTOR& kselector);

/**
    Function: WriteKmerFrequencies default implementation

    A version of <WriteKmerFrequencies()> that writes out all kmers (except unique
    ones if include_unique is false).
*/
template <class KSHAPE>
void
WriteKmerFrequencies( const vecbasevector& seqs,  const String& filename, 
                      bool include_unique = false );

/**
    Function: WriteKmerFrequencies trusted

    A version of <WriteKmerFrequencies() that writes out only those kmers
    that have at least one trusted <occurrence>, that is, an occurrence
    all of whose bases are <trusted bases>.
    It can also be configured to require at least one trusted occurrence of the kmer
    in each direction (that is, to get one trusted read of the kmer from each strand).

    Template parameters:

       KSHAPE - the <shape of kmers> to extract

    Input parameters:

       seqs - the reads
       trusted - for each base of each read, whether that base is <trusted>
       requireTrustedBothDirs - if false, just one trusted occurrence of
         a kmer is required to write out the kmer; if true, two trusted
	 occurrences -- one on each strand -- are required.
       include_unique - whether to include in the output kmers that occur
          only once

    Output parameters:

       filename - the file to which to output the frequencies of the selected
          kmers
	  
*/
template <class KSHAPE>
void
WriteKmerFrequencies( const vecbasevector& seqs,  const String& filename, 
		      const vecbitvector& trusted, 
		      int threshold = 1,
		      bool requireTrustedBothDirs = false,
		      bool include_unique = false );

/**
    Function: WriteKmerTrustedFrequencies

    A version of <WriteKmerFrequencies() that writes out only those kmers
    that have at least one trusted <occurrence>, that is, an occurrence
    all of whose bases are <trusted bases>.

    Template parameters:

       KSHAPE - the <shape of kmers> to extract

    Input parameters:

       seqs - the reads
       trusted - for each base of each read, whether that base is <trusted>

    Output parameters:

       filename - the file to which to output the frequencies of the selected
          trusted kmers
	  
*/
template <class KSHAPE>
void
WriteKmerTrustedFrequencies( const vecbasevector& seqs,  const String& filename, 
			     const vecbitvector& trusted );

#endif
// #ifndef __INCLUDE_kmer_freq_WriteKmerFrequencies_h

