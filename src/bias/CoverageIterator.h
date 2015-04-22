///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * Coverage: A simple class representing the coverage at a single
 * reference location.
 *
 * AbstractCoverageIterator: An abstract base class that is the common
 * ancestor of all the other iterators in this file.
 * 
 * SavedCoverageIterator: A class to enable serialization and deserialization
 * of other AbstractCoverageIterators.
 *
 * CoverageIterator: An iterator that returns the Coverage objects
 * implied by a reference and a sorted BAM file of mapped reads,
 * in reference-coordinate order.
 *
 * MultiCoverageIterator: An iterator that returns the Coverage objects
 * resulting from summing the coverage of multiple BAM files.
 *
 */

// Original author: Michael G. Ross <mgross@broadinstitute.org>

#ifndef BC_COVERAGE_ITERATOR_H
#define BC_COVERAGE_ITERATOR_H

#include "bias/GenomeReference.h"
#include "lookup/SAM.h"
#include "random/RNGen.h"

namespace BC
{
    /* Objects representing the coverage at a reference coordinate. */
    class Coverage
    {
    public:
        size_t chr;
        size_t loc;
        int fwd;
        int rev;
        int fwd_starts;
        int rev_starts;
        int fwd_qual;
        int rev_qual;
        int fwd_mismatches;
        int rev_mismatches;
        int fwd_deletions;
        int rev_deletions;
        int fwd_insertions;
        int rev_insertions;
        int fwd_clips;
        int rev_clips;
        int fwd_skips;
        int rev_skips;

        /* Objects representing the coverage at a reference coordinate. The
           constructor initializes each object to a meaningless value. */
        Coverage(int init_format_version = 2);
        /* Serialize the object to a binary stream. */
        void writeBinary(BinaryWriter& writer) const;
        /* Read the object from a binary stream. */
        void readBinary(BinaryReader& reader);
        /* Report the serialization size. */
        static size_t externalSizeof()
        {
            return 0; // because we don't always read the same number of bytes
        }

    private:
        int format_version;
    };

    /* Abstract base class for iterators that return a coordinate-ordered series
       of Coverages for a BAM file of aligned reads. */
    class AbstractCoverageIterator
    {
    public:
        /* Are there any more Coverages available? */
        virtual bool hasNext() = 0;
        /* Return the next Coverage in coordinate-sorted order. */
        virtual Coverage next() = 0;
        /* The number of reads the iterator has encountered so far. */
        virtual longlong numReads() = 0;
        /* The number of reads the iterator has encountered and used (or will
           use) to compute coverage so far. */
        virtual longlong numReadsUsed() = 0;
        /* Rewind the iterator to the beginning. */
        virtual void rewind() = 0;
        /* Report the libraries. */
        virtual vec<String> getLibraries() = 0;
        /* Report the production units. */
        virtual vec<String> getProductionUnits() = 0;
        /* Report the samples. */
        virtual vec<String> getSamples() = 0;
        /* Destructor. */
        virtual ~AbstractCoverageIterator();
    };

    /* A class for serializing/deserializing the stream of data produced
       by an AbstractCoverageIterator. This allows future processing to
       be much faster and more efficient than the first pass which had
       to build the coverage map from a SAM/BAM file. */
    class SavedCoverageIterator : public AbstractCoverageIterator
    {
    public:
        /* Create an empty write-only SavedCoverageIterator. */
        SavedCoverageIterator(String filename, String init_ref_fasta,
            vec<String> init_libraries, vec<String> init_production_units,
            vec<String> init_samples);
        /* Load stored coverage data. */
        SavedCoverageIterator(String filename);
        /* Destructor. */
        virtual ~SavedCoverageIterator();
        /* Add coverage data to the end of a write-only iterator. */
        void addCoverage(const Coverage& cov, longlong reads_seen,
            longlong reads_used);
        /* Are there any more Coverages available? */
        virtual bool hasNext();
        /* Return the next Coverage in coordinate-sorted order. */
        virtual Coverage next();
        /* The number of reads the iterator has encountered so far. */
        virtual longlong numReads();
        /* The number of reads the iterator has encountered and used (or will
           use) to compute coverage so far. */
        virtual longlong numReadsUsed();
        /* Rewind the iterator to the beginning. */
        virtual void rewind();
        /* Report the reference file. */
        String getReference();
        /* Report the libraries. */
        virtual vec<String> getLibraries();
        /* Report the production units. */
        virtual vec<String> getProductionUnits();
        /* Report the samples. */
        virtual vec<String> getSamples();
    private:
        String ref_fasta;
        vec<String> libraries;
        vec<String> production_units;
        vec<String> samples;
        longlong stored_reads_seen;
        longlong stored_reads_used;
        BinaryWriter* sci_writer;
        BinaryReader* sci_reader;
        const String filename;
        static const String FORMATSTRING_PREFIX;
        static const String FORMATSTRING_SUFFIX;
        int format_version;
    };

    /* Iterators that return a coordinate-ordered series of Coverages for
       a BAM file of aligned reads. */
    class CoverageIterator : public AbstractCoverageIterator
    {
    private:
        // iterator type
        const bool ignore_duplicates;
        const bool ignore_nonpf;
        const double exclude_mapq_lte;
        const double frac_to_use;
        const int pair_filtering;
        const bool only_read1;
        const bool only_read2;
        const int coverage_cap;
        const bool use_original_qscores;
        
        // genome information
        vec<String> sequence_names;
        vec<size_t> sequence_ids;
        const vecbasevector& genome;
        const vecbitvector& genome_amb;
        const vec<String>& genome_names;

        // random number generator
        const RNGen init_rand_numgen;
        RNGen rand_numgen;
        
        // position
        size_t current_chr;
        size_t current_loc;
        size_t next_chr;
        size_t next_loc;
        
        // read counts
        longlong nall_reads;
        longlong nreads;

        // read information
        bool records_remain; 
        SAM::Record next_rec;
        std::map<std::string,bool> pair_included_map;
        vec<String> libraries;
        vec<String> production_units;
        vec<String> samples;
        
        // SAM/BAM reading support
        const String bamfile;
        Logger logger;
        SAM::BAMFile* sf;
        bool secondary_alignment_seen;
       
        // coverage data
        std::deque<int> fwd_cov_queue;
        std::deque<int> rev_cov_queue;
        std::deque<int> fwd_starts_queue;
        std::deque<int> rev_starts_queue;
        std::deque<int> fwd_qual_queue;
        std::deque<int> rev_qual_queue;
        std::deque<int> fwd_mismatches_queue;
        std::deque<int> rev_mismatches_queue;
        std::deque<int> fwd_deletions_queue;
        std::deque<int> rev_deletions_queue;
        std::deque<int> fwd_insertions_queue;
        std::deque<int> rev_insertions_queue;
        std::deque<int> fwd_clips_queue;
        std::deque<int> rev_clips_queue;
        std::deque<int> fwd_skips_queue;
        std::deque<int> rev_skips_queue;
        
    public:
        /* Iterators that return a coordinate-ordered series of Coverages for
           a BAM file of aligned reads. */
        CoverageIterator(const String& bamfile,
                         const GenomeReference& gen_ref,
                         bool init_ignore_duplicates = true,
                         bool init_ignore_nonpf = true,
                         double init_exclude_mapq_lte = -1,
                         double init_frac_to_use = 1,
                         unsigned int rand_seed = 1,
                         int init_pair_filtering = 0,
                         bool init_only_read1 = false,
                         bool init_only_read2 = false,
                         int init_coverage_cap = 0,
                         bool init_use_original_qscores = true);
        /* Destructor. */
        virtual ~CoverageIterator();
        /* The number of reads the iterator has encountered so far. */
        virtual longlong numReads();
        /* The number of reads the iterator has encountered and used (or will
           use) to compute coverage so far. */
        virtual longlong numReadsUsed();
        /* Are there any more Coverages available from this BAM file? */
        virtual bool hasNext();
        /* Return the next Coverage in coordinate-sorted order. */
        virtual Coverage next();
        /* Rewind the iterator to the beginning. */
        virtual void rewind();
        /* Report the libraries. */
        virtual vec<String> getLibraries();
        /* Report the production units. */
        virtual vec<String> getProductionUnits();
        /* Report the samples. */
        virtual vec<String> getSamples();
    private:
        /* Update the internal data structure representing the coverage
           stream - ensure that the location at the top of the queue
           doesn't overlap the position of any future reads.*/
        void updateQueues();
        /* Get the next read from the BAM file. */
        bool fetchNextRead();
        /* Incorporate the next read into the coverage map. */
        void incorporateNextRead();
        /* Determine if this record should be skipped based on a number of
           settings. */
        bool skip_record(const SAM::Record& rec);
        /* Initialize the random number generator to its starting position. */
        void initRandGen();
        /* Initialize the SAM reader to its starting position. */
        void initSAM();
        /* Create a lookup for sequence names and their genome indices. */
        void map_sequence_names();
    };

    /* A coverage iterator over a set of BAM/SCI files, which returns Coverages
       representing the sum of their coverage of the genome. */
    class MultiCoverageIterator : public AbstractCoverageIterator
    {
    private:
        vec<AbstractCoverageIterator*> cov_iters;
        vec<String> sequence_names;
        vec<size_t> sequence_ids;
    public:
        /* Constructor - optional arguments only apply to BAM files in the
           list. */
        MultiCoverageIterator(const vec<String>& cov_files,
            const GenomeReference& gen_ref,
            bool ignore_duplicates = true,
            bool ignore_nonpf = true,
            double exclude_mapq_lte = -1,
            const vec<double>& frac_to_use = vec<double>(1, 1),
            bool rand_clock_seed = true,
            int pair_filtering = 0,
            bool only_read1 = false,
            bool only_read2 = false,
            int coverage_cap = 0,
            bool use_original_qscores = true);
        /* Destructor - deallocates the component CoverageIterators. */
        virtual ~MultiCoverageIterator();
        /* Are there any more Coverages available from these BAM files? */
        virtual bool hasNext();
        /* Return the next Coverage in coordinate-sorted order. */
        virtual Coverage next();
        /* The number of reads the iterator has encountered so far. */
        virtual longlong numReads();
        /* The number of reads the iterator has encountered and used (or will
           use) to compute coverage so far. */
        virtual longlong numReadsUsed();
        /* Rewind the iterator to the beginning. */
        virtual void rewind();
        /* Report the libraries. */
        virtual vec<String> getLibraries();
        /* Report the production units. */
        virtual vec<String> getProductionUnits();
        /* Report the samples. */
        virtual vec<String> getSamples();
    };

    /* A helper function for creating a coverage iterator from a BAM or SCI
       file. Note that the optional parameters only apply to BAMs. */
    AbstractCoverageIterator* load_coverage_iter(const String& cov_filename,
        const GenomeReference& gref, bool ignore_duplicates = true,
        bool ignore_nonpf = true, double exclude_mapq_lte = -1,
        double frac_to_use = 1, unsigned int rand_seed = 1,
        int pair_filtering = 0, bool only_read1 = false, 
        bool only_read2 = false, int coverage_cap = 0,
        bool use_original_qscores = true);
        
    /* Find the BAM in a Picard lane (formatted as 
       FLOWCELL.LANE[.[LIBRARY][.RUN]]). */
    String locate_lane_bam(const String& raw_lane);
    
    /* Find the BAM/SCI files, specified either via a list of BAM/SCI files
       or a list of Picard lanes (in FLOWCELL.LANE[.[LIBRARY][.RUN]] format)
       or a mixed list of the two types. */
    vec<String> locate_cov_files(String cov_param);
    
    /* Break up a file into one string per line. */
    vec<String> process_file_list(String filename);
}

SELF_SERIALIZABLE(BC::Coverage);

#endif // BC_COVERAGE_ITERATOR_H
