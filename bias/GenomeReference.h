///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * GenomeReference: A wrapper class for reading a FASTA and storing the
 * names, sequences, and ambiguous bases using quasi-appropriate AllPathsian
 * data structures.
 *
 */

// Original author: Michael G. Ross <mgross@broadinstitute.org>

#ifndef BC_GENOME_REFERENCE_H
#define BC_GENOME_REFERENCE_H

#include "Basevector.h"
#include "Bitvector.h"
#include "String.h"

namespace BC
{
    class GenomeReference
    {
    public:
        /* Load genome information from a reference FASTA. */
        GenomeReference(const String& init_ref_fasta,
            const String& cache_dir = "");
        /* Access the bases (N's are represented by a random base). */
        const vecbasevector& getBases() const;
        /* Access bitvectors indicating ambiguous bases. */
        const vecbitvector& getAmbiguities() const;
        /* Access contig names. */
        const vec<String>& getNames() const;
    private:
        vecbasevector bases;
        vecbitvector ambiguities;
        vec<String> names;
        String ref_fasta;
    };
        
    /* Look up a reference file using a Picard directory's params.txt files
       The argument is a vector of Picard-located BAMs and the function returns
       an empty string on failure. */
    String locate_reference_from_lane(const String& picard_bam);

    /* Look up a reference file using a SCI file header. */
    String locate_reference_from_sci(const String& sci_filename);
    
    /* Look up a reference file using a BAM file header. */
    String locate_reference_from_bam(const String& bam_filename);
    
    /* Look up a reference file using BAM/SCI file headers - check that all 
       refer to the same FASTA. Return empty string on failure. */
    String locate_reference(const vec<String>& cov_files);

}

#endif // BC_GENOME_REFERENCE_H
