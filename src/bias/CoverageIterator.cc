///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Original author: Michael G. Ross <mgross@broadinstitute.org>

#include "FastIfstream.h"
#include "TokenizeString.h"
#include "VecUtilities.h"
#include "bias/CoverageIterator.h"
#include "lookup/SAM2CRD.h"

namespace BC
{
    /* Objects representing the coverage at a reference coordinate. The
       constructor initializes each object to a meaningless value. */
    Coverage::Coverage(int init_format_version)
    {
        format_version = init_format_version;
        chr = numeric_limits<size_t>::max();
        loc = numeric_limits<size_t>::max();
        fwd = -1;
        rev = -1;
        fwd_starts = -1;
        rev_starts = -1;
        fwd_qual = -1;
        rev_qual = -1;
        fwd_mismatches = -1;
        rev_mismatches = -1;
        fwd_deletions = -1;
        rev_deletions = -1;
        fwd_insertions = -1;
        rev_insertions = -1;
        fwd_clips = -1;
        rev_clips = -1;
        fwd_skips = -1;
        rev_skips = -1;
        return;
    }

    /* Serialize the object to a binary stream. */
    void Coverage::writeBinary(BinaryWriter& writer) const
    {
        writer.write(chr);
        writer.write(loc);
        writer.write(fwd);
        writer.write(rev);
        writer.write(fwd_starts);
        writer.write(rev_starts);
        writer.write(fwd_qual);
        writer.write(rev_qual);
        writer.write(fwd_mismatches);
        writer.write(rev_mismatches);
        writer.write(fwd_deletions);
        writer.write(rev_deletions);
        writer.write(fwd_insertions);
        writer.write(rev_insertions);
        writer.write(fwd_clips);
        writer.write(rev_clips);
        writer.write(fwd_skips);
        writer.write(rev_skips);
    }
    
    /* Read the object from a binary stream. */
    void Coverage::readBinary(BinaryReader& reader)
    {
        reader.read(&chr);
        reader.read(&loc);
        reader.read(&fwd);
        reader.read(&rev);
        reader.read(&fwd_starts);
        reader.read(&rev_starts);
        reader.read(&fwd_qual);
        reader.read(&rev_qual);
        reader.read(&fwd_mismatches);
        reader.read(&rev_mismatches);
        reader.read(&fwd_deletions);
        reader.read(&rev_deletions);
        reader.read(&fwd_insertions);
        reader.read(&rev_insertions);
        reader.read(&fwd_clips);
        reader.read(&rev_clips);
        
        if (format_version > 1)
        {
            reader.read(&fwd_skips);
            reader.read(&rev_skips);
        }

        return;
    }

    /* Destructor. */
    AbstractCoverageIterator::~AbstractCoverageIterator()
    {
        return;
    }

    const String SavedCoverageIterator::FORMATSTRING_PREFIX =
        "::##%BC::SAVEDCOVERAGEITERATOR_FORMATv";
    const String SavedCoverageIterator::FORMATSTRING_SUFFIX = "::$$/\\)(?";

    /* Create an empty write-only SavedCoverageIterator. */
    SavedCoverageIterator::SavedCoverageIterator(String init_filename,
        String init_ref_fasta, vec<String> init_libraries,
        vec<String> init_production_units, vec<String> init_samples) :
        filename(init_filename)
    {
        ref_fasta = init_ref_fasta;
        libraries = init_libraries;
        production_units = init_production_units;
        samples = init_samples;
        stored_reads_seen = 0;
        stored_reads_used = 0;
        format_version = 2;
        
        sci_reader = NULL;
        sci_writer = new BinaryWriter(filename.c_str());
        sci_writer->write(FORMATSTRING_PREFIX + "2" + FORMATSTRING_SUFFIX);
        sci_writer->write(ref_fasta);
        sci_writer->write(samples);
        sci_writer->write(libraries);
        sci_writer->write(production_units);
        
        return;    
    }

    /* Load stored coverage data. */
    SavedCoverageIterator::SavedCoverageIterator(String init_filename) :
        filename(init_filename)
    {
        ForceAssert(IsRegularFile(filename));
        sci_reader = NULL;
        sci_writer = NULL;
        rewind();
        return;
    }

    /* Destructor. */
    SavedCoverageIterator::~SavedCoverageIterator()
    {
        if (sci_writer != NULL)
        {
            sci_writer->close();
            delete sci_writer;
        }
        
        if (sci_reader != NULL)
        {
            delete sci_reader;
        }
        
        return;
    }

    /* Add coverage data to the end of the iterator. */
    void SavedCoverageIterator::addCoverage(const Coverage& cov,
        longlong reads_seen, longlong reads_used)
    {
        if (sci_writer == NULL)
        {
            std::cerr << "Attempted to write to a read-only "
                "SavedCoverageIterator" << std::endl;
            CRD::exit(EXIT_FAILURE);
        }

        sci_writer->write(cov);
        sci_writer->write(reads_seen);
        sci_writer->write(reads_used);

        return;
    }

    /* Are there any more Coverages available? */
    bool SavedCoverageIterator::hasNext()
    {
        if (sci_reader == NULL)
        {
            std::cerr << "hasNext() is illegal when writing a "
                "SavedCoverageIterator" << std::endl;
            CRD::exit(EXIT_FAILURE);
        }
        return !sci_reader->atEOF();
    }
    
    /* Return the next Coverage in coordinate-sorted order. */
    Coverage SavedCoverageIterator::next()
    {
        if (sci_reader == NULL)
        {
            std::cerr << "next() is illegal when writing a "
                "SavedCoverageIterator" << std::endl;
            CRD::exit(EXIT_FAILURE);
        }

        if (hasNext())
        {
            Coverage next_cov(format_version);
            sci_reader->read(&next_cov);
            sci_reader->read(&stored_reads_seen);
            sci_reader->read(&stored_reads_used);
            return next_cov;
        }
        else
        {
            Coverage uninterpretable;
            return uninterpretable;
        }
    }
    
    /* The number of reads the iterator has encountered so far. */
    longlong SavedCoverageIterator::numReads()
    {
        if (sci_reader == NULL)
        {
            std::cerr << "numReads() is illegal when writing a "
                "SavedCoverageIterator" << std::endl;
            CRD::exit(EXIT_FAILURE);
        }

        return stored_reads_seen;
    }
    
    /* The number of reads the iterator has encountered and used (or will
       use) to compute coverage so far. */
    longlong SavedCoverageIterator::numReadsUsed()
    {
        if (sci_reader == NULL)
        {
            std::cerr << "numReadsUsed() is illegal when writing a "
                "SavedCoverageIterator" << std::endl;
            CRD::exit(EXIT_FAILURE);
        }

        return stored_reads_used;
    }

    /* Rewind the iterator to the beginning. */
    void SavedCoverageIterator::rewind()
    {
        if (sci_writer != NULL)
        {
            std::cerr << "rewind() is illegal when writing a "
                "SavedCoverageIterator" << std::endl;
            CRD::exit(EXIT_FAILURE);
        }

        if (sci_reader != NULL)
        {
            delete sci_reader;
        }
        
        sci_reader = new BinaryReader(filename.c_str());
        sci_reader->read(&ref_fasta);
        
        if (ref_fasta.StartsWith(FORMATSTRING_PREFIX))
        {
            format_version = ref_fasta.After(FORMATSTRING_PREFIX)
                .Before(FORMATSTRING_SUFFIX).Int();
            sci_reader->read(&ref_fasta);
            sci_reader->read(&samples);
        }
        else
        {
            format_version = 0;
            samples.push_back("<unknown>");
        }
        
        sci_reader->read(&libraries);
        sci_reader->read(&production_units);
        stored_reads_seen = 0;
        stored_reads_used = 0;

        return;
    }
    
    /* Report the reference file. */
    String SavedCoverageIterator::getReference()
    {
        return ref_fasta;
    }
    
    /* Report the libraries. */
    vec<String> SavedCoverageIterator::getLibraries()
    {
        return libraries;
    }
    
    /* Report the production units. */
    vec<String> SavedCoverageIterator::getProductionUnits()
    {
        return production_units;
    }

    /* Report the samples. */
    vec<String> SavedCoverageIterator::getSamples()
    {
        return samples;
    }
    
    /* Iterators that return a coordinate-ordered series of Coverages for
       a BAM file of aligned reads. */
    CoverageIterator::CoverageIterator(const String& init_bamfile,
                                       const GenomeReference& gen_ref,
                                       bool init_ignore_duplicates,
                                       bool init_ignore_nonpf,
                                       double init_exclude_mapq_lte,
                                       double init_frac_to_use,
                                       unsigned int init_rand_seed,
                                       int init_pair_filtering,
                                       bool init_only_read1,
                                       bool init_only_read2,
                                       int init_coverage_cap,
                                       bool init_use_original_qscores) :
        ignore_duplicates(init_ignore_duplicates),
        ignore_nonpf(init_ignore_nonpf),
        exclude_mapq_lte(init_exclude_mapq_lte),
        frac_to_use(init_frac_to_use),
        pair_filtering(init_pair_filtering),
        only_read1(init_only_read1),
        only_read2(init_only_read2),
        coverage_cap(init_coverage_cap),
        use_original_qscores(init_use_original_qscores),
        genome(gen_ref.getBases()),
        genome_amb(gen_ref.getAmbiguities()),
        genome_names(gen_ref.getNames()),
        init_rand_numgen(init_rand_seed),
        bamfile(init_bamfile),
        logger(std::cerr)
    {
        sf = NULL;
        secondary_alignment_seen = false;
        rewind();
        return;
    }
    
    /* Destructor. */
    CoverageIterator::~CoverageIterator()
    {
        delete sf;
        return;
    }
    
    /* The number of reads the iterator has encountered so far. */
    longlong CoverageIterator::numReads()
    {
        return nall_reads;
    }

    /* The number of reads the iterator has encountered and used (or will
       use) to compute coverage so far. */
    longlong CoverageIterator::numReadsUsed()
    {
        return nreads;
    }

    /* Are there any more Coverages available from this BAM file? */
    bool CoverageIterator::hasNext()
    {
        return (current_chr < genome.size() &&
                current_loc < genome[current_chr].size());
    }
    
    /* Return the next Coverage in coordinate-sorted order. */
    Coverage CoverageIterator::next()
    {
        Coverage cov; // default value is intentionally uninterpretable
        
        if (fwd_cov_queue.size() > 0)
        {
            cov.chr = current_chr;
            cov.loc = current_loc;
            cov.fwd = fwd_cov_queue.front();
            cov.rev = rev_cov_queue.front();
            cov.fwd_starts = fwd_starts_queue.front();
            cov.rev_starts = rev_starts_queue.front();
            cov.fwd_qual = fwd_qual_queue.front();
            cov.rev_qual = rev_qual_queue.front();
            cov.fwd_mismatches = fwd_mismatches_queue.front();
            cov.rev_mismatches = rev_mismatches_queue.front();
            cov.fwd_deletions = fwd_deletions_queue.front();
            cov.rev_deletions = rev_deletions_queue.front();
            cov.fwd_insertions = fwd_insertions_queue.front();
            cov.rev_insertions = rev_insertions_queue.front();
            cov.fwd_clips = fwd_clips_queue.front();
            cov.rev_clips = rev_clips_queue.front();
            cov.fwd_skips = fwd_skips_queue.front();
            cov.rev_skips = rev_skips_queue.front();
            
            fwd_cov_queue.pop_front();
            rev_cov_queue.pop_front();
            fwd_starts_queue.pop_front();
            rev_starts_queue.pop_front();
            fwd_qual_queue.pop_front();
            rev_qual_queue.pop_front();
            fwd_mismatches_queue.pop_front();
            rev_mismatches_queue.pop_front();
            fwd_deletions_queue.pop_front();
            rev_deletions_queue.pop_front();
            fwd_insertions_queue.pop_front();
            rev_insertions_queue.pop_front();
            fwd_clips_queue.pop_front();
            rev_clips_queue.pop_front();
            fwd_skips_queue.pop_front();
            rev_skips_queue.pop_front();
        }
        else if (current_chr < genome.size() &&
                 current_loc < genome[current_chr].size())
        {
            cov.chr = current_chr;
            cov.loc = current_loc;
            cov.fwd = 0;
            cov.rev = 0;
            cov.fwd_starts = 0;
            cov.rev_starts = 0;
            cov.fwd_qual = 0;
            cov.rev_qual = 0;
            cov.fwd_mismatches = 0;
            cov.rev_mismatches = 0;
            cov.fwd_deletions = 0;
            cov.rev_deletions = 0;
            cov.fwd_insertions = 0;
            cov.rev_insertions = 0;
            cov.fwd_clips = 0;
            cov.rev_clips = 0;
            cov.fwd_skips = 0;
            cov.rev_skips = 0;
        }

        if (current_loc < genome[current_chr].size() - 1)
        {
            current_loc++;
        }
        else if (current_chr < genome[current_chr].size())
        {
            // note this makes the stop-state when current_chr
            // points past the last chromosome
            current_chr++;
            current_loc = 0;
        }
        updateQueues();

        return cov;
    }
    
    /* Rewind the iterator to the beginning. */
    void CoverageIterator::rewind()
    {        
        current_chr = 0;
        current_loc = 0;
        next_chr = 0;
        next_loc = 0;

        nall_reads = 0;
        nreads = 0;

        fwd_cov_queue.clear();
        rev_cov_queue.clear();
        fwd_starts_queue.clear();
        rev_starts_queue.clear();
        fwd_qual_queue.clear();
        rev_qual_queue.clear();
        fwd_mismatches_queue.clear();
        rev_mismatches_queue.clear();
        fwd_deletions_queue.clear();
        rev_deletions_queue.clear();
        fwd_insertions_queue.clear();
        rev_insertions_queue.clear();
        fwd_clips_queue.clear();
        rev_clips_queue.clear();
        fwd_skips_queue.clear();
        rev_skips_queue.clear();

        // initialize random number generator
        initRandGen();

        // initialize SAM/BAM reading
        initSAM();

        // map the sequence names
        map_sequence_names();

        records_remain = fetchNextRead();
        updateQueues();
        return;
    }
    
    /* Report the libraries. */
    vec<String> CoverageIterator::getLibraries()
    {
        return libraries;
    }
    
    /* Report the production units. */
    vec<String> CoverageIterator::getProductionUnits()
    {
        return production_units;    
    }
    
    /* Report the samples. */
    vec<String> CoverageIterator::getSamples()
    {
        return samples;
    }
    
    /* Update the internal data structure representing the coverage
       stream - ensure that the location at the top of the queue
       doesn't overlap the position of any future reads.*/
    void CoverageIterator::updateQueues()
    {
        if (fwd_cov_queue.size() != 0 &&
            (next_chr > current_chr || next_loc > current_loc))
        {
            return; // do nothing, the next read doesn't apply yet
        }
        else
        {
            while (records_remain && next_chr == current_chr &&
                   next_loc == current_loc)
            {
                // incorporate new reads until the top of the queue no
                // longer overlaps with any read
                    
                // add coverage of next read
                incorporateNextRead();
                    
                // get a new next read
                records_remain = fetchNextRead();
            }
        }
        return;
    }

    /* Get the next read from the BAM file. */
    bool CoverageIterator::fetchNextRead()
    {
        // get the next mapped record
        while (sf->nextRecord(next_rec))
        {
            if (next_rec.isPaired())
            {
                if (next_rec.isFirstReadOfPair())
                {
                    nall_reads++;
                }
            }
            else
            {
                nall_reads++;
            }

            if (!skip_record(next_rec) &&
                (!only_read1 || next_rec.isFirstReadOfPair()) &&
                (!only_read2 || next_rec.isSecondReadOfPair()) &&
                next_rec.isMapped())
            {
                // locate the relevant chromosome
                String seqname = next_rec.getRefName();
                if (seqname.Contains(" "))
                {
                    seqname = seqname.Before(" ");
                }
                int p = BinPosition(sequence_names, seqname);
                if (p < 0)
                {
                    std::cerr << "Unable to locate contig named " << seqname
                              << " in the reference." << std::endl;
                    CRD::exit(EXIT_FAILURE);
                }
        
                next_chr = sequence_ids[p];
                next_loc = next_rec.getRefPos() - 1;

                // check that sort is maintained
                if (next_chr < current_chr ||
                    (next_loc < current_loc && next_chr == current_chr))
                {
                    std::cerr << "The source BAM file was not "
                        "properly sorted." << std::endl;
                    CRD::exit(EXIT_FAILURE);
                }

                // check that this read is actually on the contig
                if (next_loc < genome[next_chr].size())
                {
                    nreads++;        
                    return true;
                }
                else
                {
                    std::cerr << "Warning: " << next_rec.getQueryName()
                        << " is being ignored because its alignment begins "
                        "after contig " << seqname << " ends\n";
                }
            }
        }

        // no valid read -> ran out of records
        return false;
    }

    /* Incorporate the next read into the coverage map. */
    void CoverageIterator::incorporateNextRead()
    {
        SAM::Alignment aln(next_rec, logger);
        size_t p1 = aln.getReadStart() - aln.getLeftHardClip();
        size_t p2 = aln.getRefStart();
        size_t required_length = aln.getRefEnd() - current_loc;
        bool rc1 = next_rec.isReversed();
        std::vector<uchar> const& qual =
            (use_original_qscores &&
            next_rec.getOriginalQualityScores().size() > 0) ?
            next_rec.getOriginalQualityScores() :
            next_rec.getQualityScores();
        std::string const& seq = next_rec.getSequence();
        
        if (required_length > fwd_cov_queue.size())
        {
            fwd_cov_queue.resize(required_length, 0);
            rev_cov_queue.resize(required_length, 0);
            fwd_starts_queue.resize(required_length, 0);
            rev_starts_queue.resize(required_length, 0);
            fwd_qual_queue.resize(required_length, 0);
            rev_qual_queue.resize(required_length, 0);
            fwd_mismatches_queue.resize(required_length, 0);
            rev_mismatches_queue.resize(required_length, 0);
            fwd_deletions_queue.resize(required_length, 0);
            rev_deletions_queue.resize(required_length, 0);
            fwd_insertions_queue.resize(required_length, 0);
            rev_insertions_queue.resize(required_length, 0);
            fwd_clips_queue.resize(required_length, 0);
            rev_clips_queue.resize(required_length, 0);
            fwd_skips_queue.resize(required_length, 0);
            rev_skips_queue.resize(required_length, 0);
        }

        // if this read is the second read of a mapped pair,
        // we shouldn't count its start point because it will
        // be correlated with its mate's start point
        bool marked_start = !(!next_rec.isPaired() ||
                              next_rec.isFirstReadOfPair() ||
                              !next_rec.isMateMapped());

        // p1 > 0 indicates the beginning of the read has been clipped
        if (p2 < genome[current_chr].size())
        {
            if (rc1)
            {
                rev_clips_queue[p2 - current_loc] += p1;
            }
            else
            {
                fwd_clips_queue[p2 - current_loc] += p1;
            }
        }
        
        for (size_t sam_block = 0; sam_block < aln.size(); sam_block++)
        {
            p1 = aln[sam_block].getReadStart() - aln.getLeftHardClip();
            p2 = aln[sam_block].getRefStart();
            if (aln[sam_block].isRefEater() && !aln[sam_block].isReadEater())
            {
                char cigar_op = aln[sam_block].getCigarOp();
                for (size_t d = p2; d < (p2 + aln[sam_block].getLength()); d++)
                {
                    if (d < genome[current_chr].size())
                    {
                        if (rc1)
                        {
                            if (cigar_op == 'D')
                            {
                                rev_deletions_queue[d - current_loc]++;
                            }
                            else if (cigar_op == 'N')
                            {
                                rev_skips_queue[d - current_loc]++;
                            }
                            else
                            {
                                FatalErr("Unexpected CIGAR op " + cigar_op);
                            }
                        }
                        else
                        {
                            if (cigar_op == 'D')
                            {
                                fwd_deletions_queue[d - current_loc]++;
                            }
                            else if (cigar_op == 'N')
                            {
                                fwd_skips_queue[d - current_loc]++;
                            }
                            else
                            {
                                FatalErr("Unexpected CIGAR op " + cigar_op);
                            }
                        }
                    }
                }
            }
            else if (aln[sam_block].isReadEater() && 
                !aln[sam_block].isRefEater())
            {
                // because insertions do not, by definition, correspond
                // precisely to a reference base, the errors are all charged
                // to the current location on the reference - so the insertion
                // can be considered to be occurring to the left of it
                if (p2 < genome[current_chr].size())
                {
                    if (rc1)
                    {
                        rev_insertions_queue[p2 - current_loc] +=
                            aln[sam_block].getLength();
                    }
                    else
                    {
                        fwd_insertions_queue[p2 - current_loc] +=
                            aln[sam_block].getLength();
                    }
                }
            }
            else if (aln[sam_block].isRefEater() &&
                aln[sam_block].isReadEater())
            {
                for (size_t x = 0; x < aln[sam_block].getLength(); x++)
                {
                    if (p2 < genome[current_chr].size())
                    {
                        if (!marked_start)
                        {
                            if (rc1)
                            {
                                rev_starts_queue[p2 - current_loc]++;
                            }
                            else
                            {
                                fwd_starts_queue[p2 - current_loc]++;
                            }
                            marked_start = true;
                        }
    
                        if ((coverage_cap == 0)
                            || (rev_cov_queue[p2 - current_loc]
                                + fwd_cov_queue[p2 - current_loc]
                                < coverage_cap))
                        {
                            if (!genome_amb[current_chr][p2])
                            {
                                if (Base::val2Char(genome[current_chr][p2])
                                    != seq[p1])
                                {
                                    if (rc1)
                                    {
                                        rev_mismatches_queue[p2 - 
                                            current_loc]++;
                                    }
                                    else
                                    {
                                        fwd_mismatches_queue[p2 - 
                                            current_loc]++;
                                    }
                                }
                            }
                            if (rc1)
                            {
                                rev_cov_queue[p2 - current_loc]++;
                            }
                            else
                            {
                                fwd_cov_queue[p2 - current_loc]++;
                            }
                            
                            // if some reads over this location lacked qscores,
                            // flag the location and don't record a qscore sum
                            if (qual.size() == 0)
                            {
                                if (rc1)
                                {
                                    rev_qual_queue[p2 - current_loc] = -1;
                                }
                                else
                                {
                                    fwd_qual_queue[p2 - current_loc] = -1;
                                }
                            }
                            else
                            {
                                if (rc1 && rev_qual_queue[p2 - current_loc] != 
                                    -1)
                                {
                                    rev_qual_queue[p2 - current_loc] += 
                                        qual[p1];
                                }
                                else if (fwd_qual_queue[p2 - current_loc] != -1)
                                {
                                    fwd_qual_queue[p2 - current_loc] += 
                                        qual[p1];
                                }
                            }
                        }
                    }
                    ++p1;
                    ++p2;
                }
            }
        }

        // Pos1() < seq.size() indicates that the end of the read was clipped,
        // we charge this to the last base that the read mapped to in order
        // to keep within the current bounds of the queues (and to avoid
        // complicating other parts of the program by extending them in this
        // instance)
        if (p2 >= 1 && p2 < genome[current_chr].size() + 1)
        {
            if (rc1)
            {
                rev_clips_queue[p2 - 1 - current_loc] +=
                    (seq.size() - aln.getReadEnd() + aln.getLeftHardClip());
            }
            else
            {
                fwd_clips_queue[p2 - 1 - current_loc] +=
                    (seq.size() - aln.getReadEnd() + aln.getLeftHardClip());
            }
        }
        
        return;        
    }

    /* Determine if this record should be skipped based on a number of
       settings. */
    bool CoverageIterator::skip_record(const SAM::Record& rec)
    {
        if (!rec.isPrimary() && !secondary_alignment_seen)
        {
            std::cerr << "CoverageIterator warning: " << bamfile
                << " contains secondary alignments, which are ignored."
                << std::endl;
            secondary_alignment_seen = true;
        }
        
        if (!rec.isPrimary() ||
            !rec.isMapped() ||
            (ignore_duplicates && rec.isDuplicate()) ||
            (ignore_nonpf && !rec.isPF()) ||
            (rec.getMapQ() <= exclude_mapq_lte))
        {    
            return true;
        }
    
        // if pair filtering is enabled, ignore read pairs not mapped to
        // the same chromosome, or with inserts greater than pair_filtering
        if (rec.isPaired() && pair_filtering > 0)
        {
            if ((rec.isMapped() != rec.isMateMapped()) ||
                (rec.getRefName() != rec.getMateRefName()) ||
                (abs(rec.getInsertSize()) > pair_filtering))
            {
                return true;
            }
        }

        // random downsampling    
        if (frac_to_use < 1)
        {
            if (rec.isPaired()) // behave consistently with mates
            {
                // check if we've already included/discarded
                // another segment from this template, if so, treat
                // this read identically
                std::string qname = rec.getQueryName();
                std::map<std::string,bool>::iterator pi_iter =
                    pair_included_map.find(qname);
                if (pi_iter != pair_included_map.end())
                {
                    // act as your mate did
                    return pi_iter->second;
                }
                else // mate not previously encountered, save choice
                {
                    bool skip_read = rand_numgen.next() >
                        frac_to_use * RNGen::RNGEN_RAND_MAX;
                    pair_included_map.insert
                        (std::pair<std::string,bool>(qname, skip_read));
                    return skip_read;
                }
            }
            else // no mate exists, act randomly
            {
                return (rand_numgen.next() >
                    frac_to_use * RNGen::RNGEN_RAND_MAX);
            }
        }

        return false;
    }

    /* Initialize the random number generator to its starting position. */
    void CoverageIterator::initRandGen()
    {
        rand_numgen = init_rand_numgen;
        pair_included_map.clear();
        return;
    }
    
    /* Initialize the SAM reader to its starting position. */
    void CoverageIterator::initSAM()
    {
        if (sf != NULL)
        {
            delete sf;
        }
                
        sf = new SAM::BAMFile(bamfile);
        
        const std::vector<SAM::ReadGroup>& read_groups = sf->getReadGroups();
        for (size_t r = 0; r < read_groups.size(); r++)
        {
            if (read_groups[r].getLibrary().empty())
            {
                libraries.push_back("<unknown>");
            }
            else
            {
                libraries.push_back(read_groups[r].getLibrary());
            }
            
            if (read_groups[r].getPlatformUnit().empty())
            {
                production_units.push_back("<unknown>");
            }
            else
            {
                production_units.push_back(read_groups[r].getPlatformUnit());
            }
            
            if (read_groups[r].getSample().empty())
            {
                samples.push_back("<unknown>");
            }
            else
            {
                samples.push_back(read_groups[r].getSample());
            }
        }
        UniqueSort(libraries);
        UniqueSort(production_units);
        UniqueSort(samples);

        return;
    }

    /* Create a lookup for sequence names and their genome indices. */
    void CoverageIterator::map_sequence_names()
    {
        if ( !validateReferenceDictionary(*sf,genome_names,genome) )
        {
            std::cerr << "Genome data and SAM header don't agree.\n";
            CRD::exit(EXIT_FAILURE);
        }
        sequence_names = genome_names;
        size_t nNames = genome_names.size();
        sequence_ids.clear();
        sequence_ids.reserve(nNames);
        for ( size_t idx = 0; idx != nNames; ++idx )
            sequence_ids.push_back(idx);
        SortSync(sequence_names, sequence_ids);
        return;
    }

    /* A coverage iterator over a set of BAM/SCI files, which returns Coverages
       representing the sum of their coverage of the genome. Optional
       arguments only apply to the BAM files. */
    MultiCoverageIterator::MultiCoverageIterator(const vec<String>& cov_files,
        const GenomeReference& gen_ref,
        bool ignore_duplicates,
        bool ignore_nonpf,
        double exclude_mapq_lte,
        const vec<double>& frac_to_use,
        bool rand_clock_seed,
        int pair_filtering,
        bool only_read1,
        bool only_read2,
        int coverage_cap,
        bool use_original_qscores)
    {
        if (cov_files.size() < 1)
        {
            std::cerr << "MultiCoverageIterator needs at least one input."
                      << std::endl;
            CRD::exit(EXIT_FAILURE);
        }

        cov_iters.resize(cov_files.size());
        
        if (frac_to_use.size() != 1 && frac_to_use.size() != cov_files.size())
        {
            std::cerr << "frac_to_use's length must be 1 or equal to "
                << "cov_files' length\n";
            CRD::exit(EXIT_FAILURE);
        }
        
        for (size_t b = 0; b < cov_files.size(); b++)
        {
            double frac = frac_to_use[0];
            if (frac_to_use.size() > 1)
            {
                frac = frac_to_use[b];
            }
            
            unsigned int rand_seed = b;
            if (rand_clock_seed)
            {
                timeval tval;
                gettimeofday(&tval, NULL);
                rand_seed = tval.tv_usec;
            }
                        
            cov_iters[b] = load_coverage_iter(cov_files[b], gen_ref,
                ignore_duplicates, ignore_nonpf, exclude_mapq_lte,
                frac, rand_seed, pair_filtering, only_read1,
                only_read2, coverage_cap, use_original_qscores);
        }

        return;
    }

    /* Destructor - deallocates the component CoverageIterators. */
    MultiCoverageIterator::~MultiCoverageIterator()
    {
        for (size_t i = 0; i < cov_iters.size(); i++)
        {
            delete cov_iters[i];
        }
    }

    /* Are there any more Coverages available from these BAM files? */
    bool MultiCoverageIterator::hasNext()
    {
        // it should be impossible for this to get out-of-sync
        return cov_iters[0]->hasNext();
    }

    /* Return the next Coverage in coordinate-sorted order. */
    Coverage MultiCoverageIterator::next()
    {
        Coverage combo_next = cov_iters[0]->next();
        for (size_t i = 1; i < cov_iters.size(); i++)
        {
            Coverage added_cov = cov_iters[i]->next();
            if (added_cov.chr != combo_next.chr ||
                added_cov.loc != combo_next.loc)
            {
                std::cerr << "CoverageIterators are unsynced!" << std::endl;
                CRD::exit(EXIT_FAILURE);
            }
            combo_next.fwd += added_cov.fwd;
            combo_next.rev += added_cov.rev;
            combo_next.fwd_starts += added_cov.fwd_starts;
            combo_next.rev_starts += added_cov.rev_starts;
            combo_next.fwd_qual += added_cov.fwd_qual;
            combo_next.rev_qual += added_cov.rev_qual;
            combo_next.fwd_mismatches += added_cov.fwd_mismatches;
            combo_next.rev_mismatches += added_cov.rev_mismatches;
            combo_next.fwd_deletions += added_cov.fwd_deletions;
            combo_next.rev_deletions += added_cov.rev_deletions;
            combo_next.fwd_insertions += added_cov.fwd_insertions;
            combo_next.rev_insertions += added_cov.rev_insertions;
            combo_next.fwd_clips += added_cov.fwd_clips;
            combo_next.rev_clips += added_cov.rev_clips;
            combo_next.fwd_skips += added_cov.fwd_skips;
            combo_next.rev_skips += added_cov.rev_skips;
        }
        return combo_next;
    }

    /* The number of reads the iterator has encountered so far. */
    longlong MultiCoverageIterator::numReads()
    {
        longlong nall_reads = 0;
        for (size_t i = 0; i < cov_iters.size(); i++)
        {
            nall_reads += cov_iters[i]->numReads();
        }
        return nall_reads;
    }

    /* The number of reads the iterator has encountered and used (or will
       use) to compute coverage so far. */
    longlong MultiCoverageIterator::numReadsUsed()
    {
        longlong nreads = 0;
        for (size_t i = 0; i < cov_iters.size(); i++)
        {
            nreads += cov_iters[i]->numReadsUsed();
        }
        return nreads;
    }
    
    /* Rewind the iterator to the beginning. */
    void MultiCoverageIterator::rewind()
    {
        for (size_t c = 0; c < cov_iters.size(); c++)
        {
            cov_iters[c]->rewind();
        }
        
        return;
    }

    /* Report the libraries. */
    vec<String> MultiCoverageIterator::getLibraries()
    {
        vec<String> libraries;
        for (size_t i = 0; i < cov_iters.size(); i++)
        {
            vec<String> ci_libraries = cov_iters[i]->getLibraries();
            libraries.insert(libraries.end(), ci_libraries.begin(),
                ci_libraries.end());
        }
        UniqueSort(libraries);
        return libraries;
    }
    
    /* Report the production units. */
    vec<String> MultiCoverageIterator::getProductionUnits()
    {
        vec<String> production_units;
        for (size_t i = 0; i < cov_iters.size(); i++)
        {
            vec<String> ci_production_units =
                cov_iters[i]->getProductionUnits();
            production_units.insert(production_units.end(),
                ci_production_units.begin(), ci_production_units.end());
        }
        UniqueSort(production_units);
        return production_units;
    }
    
    /* Report the samples. */
    vec<String> MultiCoverageIterator::getSamples()
    {
        vec<String> samples;
        for (size_t i = 0; i < cov_iters.size(); i++)
        {
            vec<String> ci_samples = cov_iters[i]->getSamples();
            samples.insert(samples.end(), ci_samples.begin(), ci_samples.end());
        }
        UniqueSort(samples);
        return samples;
    }

    /* A helper function for creating a coverage iterator from a BAM or SCI
       file. Note that the optional parameters only apply to BAMs. */
    AbstractCoverageIterator* load_coverage_iter(const String& cov_filename,
        const GenomeReference& gref,  bool ignore_duplicates,
        bool ignore_nonpf, double exclude_mapq_lte,
        double frac_to_use, unsigned int rand_seed,
        int pair_filtering, bool only_read1, 
        bool only_read2, int coverage_cap, bool use_original_qscores)
    {
        if (cov_filename.EndsWith(".bam") || cov_filename.EndsWith(".sam"))
        {
            return new BC::CoverageIterator(cov_filename, gref,  
                ignore_duplicates, ignore_nonpf, exclude_mapq_lte,
                frac_to_use, rand_seed, pair_filtering, only_read1, 
                only_read2, coverage_cap, use_original_qscores);
        }
        else if (cov_filename.EndsWith(".sci"))
        {
            if (frac_to_use != 1)
            {
                FatalErr("frac_to_use must be 1 for an SCI file");
            }
            return new BC::SavedCoverageIterator(cov_filename);
        }
        else
        {
            std::cerr << "Failed to load " << cov_filename << std::endl;
            CRD::exit(1);
        }
    }
    
    /* Find the BAM in a Picard lane (formatted as 
       FLOWCELL.LANE[.[LIBRARY][.RUN]]). */
    String locate_lane_bam(const String& raw_lane)
    {   
        String fc = raw_lane.Before(".");
        String lane = raw_lane.After(".");
        
        // allow for extended FLOWCELL.LANE.LIBRARY specifier,
        // useful for indexed runs that produce multiple library subdirs
        String libname = "";
        if (lane.Contains("."))
        {
            libname = lane.After(".");
            lane = lane.Before(".");
        }
        String runname = "";
        if (libname.Contains("."))
        {
            runname = libname.After(".");
            libname = libname.Before(".");   
        }
        
        String dir = "/seq/picard/" + fc;
        
        vec<String> all = AllFiles(dir);
        vec<String> run_dirs;
        for (size_t j = 0; j < all.size(); j++)
        {
            if (IsDirectory(dir + "/" + all[j]))
            {
                run_dirs.push_back(all[j]);
            }
        }
        if (!runname.empty())
        {
            if (!Member(run_dirs, runname))
            {
                std::cerr << "Can't find RUN." << std::endl;
                CRD::exit(EXIT_FAILURE);
            }
            run_dirs.clear();
            run_dirs.push_back(runname);
        }
        if (!run_dirs.solo())
        {
            std::cerr << "Found " << run_dirs.size()
                << " run directories inside "
                << "flowcell " << dir << " directory." << std::endl;
            Sort(run_dirs);
            for (size_t j = 0; j < run_dirs.size(); j++)
            {
                std::cerr << "\t" << run_dirs[j] << std::endl;
            }
            CRD::exit(EXIT_FAILURE);
        }
        dir += "/" + run_dirs.back() + "/" + lane;
        if (!IsDirectory(dir))
        {
            std::cerr << "Can't find lane subdirectory " << dir
                << std::endl;
            CRD::exit(EXIT_FAILURE);
        }
        all = AllFiles(dir);
        if (!libname.empty())
        {
            size_t lib = 0;
            String libsubdir = "";
            while (lib < all.size() && libsubdir.empty())
            {
                if (all[lib] == libname)
                {
                    libsubdir = all[lib];
                }
                lib++;
            }
            if (libsubdir.empty())
            {
                std::cerr << libname << " not found in " << dir << "."
                    << std::endl;
                CRD::exit(EXIT_FAILURE);
            }
            else
            {
                dir += "/" + libsubdir;
            }
        }
        else if (!all.solo())
        {
            std::cerr << "Found " << all.size() << " library directories."
                 << std::endl;
            std::cerr << "Please specify one of: " << std::endl;
            for (size_t lib = 0; lib < all.size(); lib++)
            {
                if (IsDirectory(dir + "/" + all[lib]))
                {
                    std::cerr << fc << "." << lane << "." << all[lib]
                        << std::endl;
                }
            }
            CRD::exit(EXIT_FAILURE);
        }
        else
        {
            dir += "/" + all[0];
        }
        
        if (!IsRegularFile(dir + "/finished.txt"))
        {
            std::cerr << "Picard directory " << dir <<
                " does not appear to be completely processed.\n";
            CRD::exit(EXIT_FAILURE);
        }
        
        String bam_filename = dir + "/" + fc + "." + lane
            + ".aligned.duplicates_marked.bam";
    
        return bam_filename;
    }
    
    /* Find the BAM/SCI files, specified either via a list of BAM/SCI files
       or a list of Picard lanes (in FLOWCELL.LANE[.[LIBRARY][.RUN]] format)
       or a mixed list of the two types. */
    vec<String> locate_cov_files(String cov_param)
    {
        bool bam_failure = false;
    
        vec<String> covs;
        if (cov_param.StartsWith("@"))
        {
            covs = process_file_list(cov_param.After("@"));
            String cov_param_dir = "";
            if (cov_param.Contains("/"))
            {
                cov_param_dir = cov_param.After("@").RevBefore("/") + "/";
            }
            for (size_t i = 0; i < covs.size(); i++)
            {
                if (covs[i].EndsWith(".bam") || covs[i].EndsWith(".sam") ||
                    covs[i].EndsWith(".sci"))
                {
                    if (!covs[i].StartsWith("/"))
                    {
                        covs[i] = cov_param_dir + covs[i];
                    }
                }
            }
        }
        else
        {
            Tokenize(cov_param, covs);
        }
        
        for (size_t i = 0; i < covs.size(); i++)
        {
            if (!covs[i].EndsWith(".bam") && !covs[i].EndsWith(".sam") &&
                !covs[i].EndsWith(".sci"))
            {
                covs[i] = locate_lane_bam(covs[i]);
            }
            if (!IsRegularFile(covs[i]))
            {
                std::cerr << "Can't find " << covs[i] << std::endl;
                CRD::exit(EXIT_FAILURE);
            }
            covs[i] = RealPath(covs[i]);
        }
    
        return covs;
    }
    
    /* Break up a file into one string per line. */
    vec<String> process_file_list(String filename)
    {
        vec<String> toks;
        String line;
        fast_ifstream in(filename);
        while (!in.fail())
        {
            getline(in, line);
            DeleteLeadingWhiteSpace(line);
            DeleteTrailingWhiteSpace(line);
    
            if (line != "")
            {
                toks.push_back(line);
            }
        }
        return toks;
    }
}
