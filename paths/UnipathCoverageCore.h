/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// UnipathCoverageCore - see UnipathCoverage.cc for description of what this
// is about.

#include "Basevector.h"
#include "CoreTools.h"
#include "paths/KmerPath.h"
#include "paths/PdfEntry.h"

// Given a unipath I, estimate its number of copies in the genome, by returning
// probabilities for each possible copy number (stopping at very small values).
// This is actually not mathematically well-defined, since it depends on some
// assumption regarding the distribution of copy numbers amongst all unipaths in
// the genome.  What we do here is treat all copy numbers "within the range of 
// interest" as occurring with equal frequency.

void CopyNumber( 
     /* inputs: */ int n,                // number of reads sharing a kmer with I
                   int Ilen,             // length of unipath, in kmers
                   double rlen,          // average read length, in kmers
                   longlong allreads,    // total number of reads
                   longlong gsize,       // genome size
                   int K,                // kmer length -- used only with errors
     /* output: */ PdfEntryVec& copyno_prob,
     /* params: */ double thresh = 0.01,
                   double error_rate = 0.0, // see "Errors" below
                   double* q_ptr = NULL );

// If unipath I has length Ilen, then the number of positions on the genome 
// at which a KmerPath overlapping I and of length rlen can start is 
// Ilen + rlen - 1.  Let q = (allreads/gsize) * ( Ilen + rlen - 1 ).  
// This is the expected number of reads that would share a kmer with the 
// unipath, if reads were uniformly disributed.  If reads are of differing 
// lengths, let rlen be the average read length; then q is still the expected 
// number of hits.
//
// Then if the copy number of I is C, and x = qC, then the probability
// that n reads share a kmer with I is
// p(n) = e-x x^n / n!, as in the Poisson pdf
//      = e-(qC) (qC)^n / n!.
//
// Therefore what CopyNumber does is to evaluate this expression for successive
// values of C, but only returns values which are at least 1% (thresh) of the 
// largest value.  These values are then normalized so that their sum is one.
//
// Errors: there may be errors in the data, in which case we would
// like to model copy number zero for kmers which are not actually in
// the genome.  Here is a crude model: We are less likely to see an
// error kmer than a true one, so we weight the prior of copy number zero
// at K*error_rate of other copy numbers.  Having seen an error, it's also
// unlikely to happen again, since that demands another error that matches
// perfectly (probability error_rate/3).  These assume that there's a
// single place in the genome where the matching errors come from, which
// is totally bogus: kmers one base off from high-frequency stuff are much
// more of a problem.  But I can't see how to take care of that.

void UnipathCoverageCore( 
     // inputs:
     const int K, const vecKmerPath& paths, const vecKmerPath& paths_rc, 
     const vec<tagged_rpint>& pathsdb, const vecKmerPath& unipaths, 
     const vec<tagged_rpint>& unipathsdb, const vec<int>& read_lengths, 
     // output:
     VecPdfEntryVec& p,
     // optional args:
     const int UNIPATH_TO_TRACE = -1, const double THRESH = 0.0001, 
     const double ERROR_RATE = 0.0, const unsigned int USE_THIS_GENOME_SIZE = 0, 
     const int PLOIDY = 1, const vec<double>* p_unipath_bias = 0 );
