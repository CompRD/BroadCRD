/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/*    Given a <.fastb> file of reads, and the table T of strong (trusted)  k-mers in */
/*    them (together with their multiplicities, as produced by */
/*    <FindKmerFrequencies>), edit the reads as follows. */
   
/*    Scan each read R, starting at the left, for strong kmers. Mark base positions that */
/*    belong to strong k-mers as strong. If no strong k-mers can be found within the original  */
/*    read then try to uncover them by mutating the bases in a following way. Scan the read from  */
/*    left to right and try 3k possible mutations for each k-mer. If only one such mutation */
/*    produces a strong k-mer then mark the corresponding base positions as strong.  */

/*    If the original or mutated read contains strong k-mers then take the leftmost one  */
/*    (under the assumption that the beginning of the read has a higher quality) and  */
/*    try to extend that kmer by including (in turn) the neighboring positions to the left  */
/*    and right. If the original base at those positions does not form an adjacent strong kmer */
/*    and was not previously marked as strong then mutate it. If only one of the 3 possible  */
/*    mutations forms the strong kmer then accept it and move to the next position. Otherwise,  */
/*    stop and move to the next seed to the right and follow the same procedure. */

/*    If no read trimming is allowed then accept the modified read if it is composed only of  */
/*    the strong k-mers and the number of mutations is lower than MAX_ERRORS as specified  */
/*    on the command line.  */

/*    If trimming is allowed, then accept an antirely strong sub-read with the highest sum of match */
/*    scores (the single base match score is 1 if the base was not mutated and -1 if it was). */
/*    In that case MAX_ERRORS is disregarded and sub-reads may contain many mutations. */

#ifndef KMERFREQ_MAKEREADCORRECTIONS_H
#define KMERFREQ_MAKEREADCORRECTIONS_H

#include "Basevector.h"
#include "Qualvector.h"
#include "kmer_freq/KmerFrequencyTable.h"
#include "paths/BaseErrorProb.h"

void MakeReadCorrections(vecbasevector& reads, 
                         const vec<int>& i_reads_todo,
                         const vecqualvector& quals,
                         const KmerShortMap& KmerMap,
                         vec<Bool>& read_is_strong,
                         const bool verbose = false,
                         const int max_errors = 1, 
                         const bool trim_reads = false,
                         const bool thorough = false,
			 const bool keep_partial = false);

// Convenience method when you want to call the above with just one value of K.
template <int K>
void MakeReadCorrections(vecbasevector& reads, 
                         const vec<int>& i_reads_todo,
                         const vecqualvector& quals,
                         vec<Bool>& read_is_strong,
                         const String& filename,
                         const bool verbose = false,
                         const int max_errors = 1,
                         const bool trim_reads = false,
                         const bool thorough = false,
			 const bool keep_partial = false);


#endif
