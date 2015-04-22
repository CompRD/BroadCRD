///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// BadCoverage. This program takes as input aligned reads, given by a BAM file.
// For each position on the genome, it either computes the sum of the qualities
// of all reads mapped to that location (the default) *or* it computes the
// "worst-strand quality coverage": the minimum over the two strands, of the sum
// of quality scores of read bases covering it (the original method).  It then
// reports the "bad fraction": the fraction of reference bases whose
// quality coverage is at least ten-fold below the mean.
//
// The "worst-strand" mode will thus identify not only places where coverage is
// low, but also places where coverage is poor.  For example it will identify
// loci where reads in one direction fail after passing through a given
// sequence.  This is a known problem for assembly.  It breaks contigs.
// However, the default mode, which is insensitive to this problem, is
// generally preferred because it produces a much more reliable result when
// relatively few (less than 1 million) reads are available, and is easier
// to analyze due to its linear behavior.
//
// -----------------------------------------------------------------------------
// SPECIAL SYSTEM REQUIREMENTS
// -----------------------------------------------------------------------------
//
// samtools needs to be in your path.
//
// -----------------------------------------------------------------------------
// REFERENCE SEQUENCE
// -----------------------------------------------------------------------------
//
// See option REF_FASTA.
//
// If no fasta file for the reference is provided as input, the program will try
// to infer the fasta file from the BAM file.  Unfortunately, BAM files do not
// in general provide a definitive description of the reference sequence.
//
// -----------------------------------------------------------------------------
// MOTIFS
// -----------------------------------------------------------------------------
//
// See options TEST_MOTIFS.
//
// A motif consists of a sequence of A, C, G, T, or N, possibly intermixed with
// "fixed length base content restrictions", for example (G+C >= 0.8)^100, which
// specifies a sequence of length 100 having GC content at least 80%.
//
// The motif syntax allows abbreviation using single levels of parentheses and
// exponents.  White space is ignored.  The following are examples of valid
// motifs:
//
// (GC)^3
// N^100 (GC)^3 N^100
// C^5 G C^4
// (G+C >= 0.8)^100
// (G >= 0.9 & G+C = 1)^10
// (G >= 0.9 | G+C = 1)^10
//
// Ambiguous bases never match a motif.

#include "Basevector.h"
#include "Bitvector.h"
#include "FastIfstream.h"
#include "FeudalMimic.h"
#include "Intvector.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "VecUtilities.h"
#include "bias/CoverageIterator.h"
#include "bias/BadMotifCore.h"
#include "bias/GenomeInfoCache.h"
#include "bias/GenomeReference.h"
#include "math/Functions.h"

namespace
{
    /* Print a coverage string for a location - N=[1-9] indicates N*10%
       coverage, ' ' indicates 0%, '+' indicates 100% or greater. */
    void print_coverage_string(vec<longlong> coverage_counts, double mean_cov,
                               unsigned int whitespace_buffer, ostream& out)
    {
        if (std::isnan(mean_cov))
        {
            return;
        }
        for (unsigned int k = 0; k < whitespace_buffer; k++)
        {
            out << ' ';
        }
        for (unsigned int k = 0; k < coverage_counts.size(); k++)
        {
            if (mean_cov == 0)
            {
                out << 'X';
            }
            else
            {
                int local_relcount = (coverage_counts[k] / mean_cov * 10 + 0.5);
                if (local_relcount > 9)
                {
                    out << '+';
                }
                else if (local_relcount == 0)
                {
                    out << '_';
                }
                else
                {
                    out << local_relcount;
                }
            }
        }
        for (unsigned int k = 0; k < whitespace_buffer; k++)
        {
            out << ' ';
        }
        return;
    }
    
    /* Print a quality string for a location - N=[1-9] indicates N*10% of mean 
        quality, ' ' indicates 0%, '+' indicates 100% or greater. */
    void print_quality_string(vec<longlong> quality_counts,
        vec<longlong> coverage_counts, double mean_qual,
        unsigned int whitespace_buffer, ostream& out)
    {
        if (std::isnan(mean_qual))
        {
            return;
        }
        for (unsigned int k = 0; k < whitespace_buffer; k++)
        {
            out << ' ';
        }
        for (unsigned int k = 0; k < coverage_counts.size(); k++)
        {
            if (coverage_counts[k] == 0)
            {
                out << 'X';
            }
            else
            {
                int local_relqual = double(quality_counts[k])
                    / double(coverage_counts[k]) / mean_qual * 10 + 0.5;
                if (local_relqual > 9)
                {
                    out << '+';
                }
                else if (local_relqual == 0)
                {
                    out << '_';
                }
                else
                {
                    out << local_relqual;
                }
            }
        }
        for (unsigned int k = 0; k < whitespace_buffer; k++)
        {
            out << ' ';
        }
        return;
    }
    
    /* Print the reference in this location - aligned with output of
       print_coverage_string. Capitals indicate the actual location,
       lowercase indicates immediate surroundings. */
    void print_genome_string(const basevector& genome, int start_index,
                             int end_index, unsigned int lowercase_buffer,
                             ostream& out)
    {
        int primer_start = start_index - lowercase_buffer;
        int primer_end = end_index + lowercase_buffer;
        for (int k = primer_start; k <= primer_end; k++)
        {
            if (k < 0 || static_cast<unsigned>(k) >= genome.size())
            {
                out << '?';
            }
            else
            {
                if (k < start_index || k > end_index)
                {
                    out << char(tolower(as_base(genome[k])));
                }
                else
                {
                    out << char(toupper(as_base(genome[k])));
                }
            }
        }
        return;
    }
    
    /* Return a vector of motif strings. */
    vec<String> parse_motifs(String test_motifs)
    {
        vec<String> motifs;
        String line;
        if (test_motifs.StartsWith("@"))
        {
            motifs = BC::process_file_list(test_motifs.After("@"));
        }
        else
        {
            ParseStringSet(test_motifs, motifs);
        }
        for (size_t m = 0; m < motifs.size(); m++)
        {
            motifs[m] = ToUpper(motifs[m]);
            DeleteLeadingWhiteSpace(motifs[m]);
            DeleteTrailingWhiteSpace(motifs[m]);
        }
        return motifs;
    }
    
    /* Load intervals (typically the "bad promoters" list). */
    void load_intervals(vec<size_t>& chr, vec<size_t>& start,
        vec<size_t>& stop, const BC::GenomeReference& gen_ref,
        String intervals)
    {
        const vec<String>& names = gen_ref.getNames();
        const vecbasevector& bases = gen_ref.getBases();
        vec<size_t> init_chr;
        vec<size_t> init_start;
        vec<size_t> init_stop;
        
        if (!IsRegularFile(intervals))
        {
            FatalErr(intervals  + " does not appear to be a regular file.");
        }
        
        fast_ifstream testin(intervals);
        String line;
        bool baitlist_format = intervals.EndsWith(".interval_list");
        while (!testin.fail())
        {
            getline(testin, line);
            if (!testin.fail())
            {
                if (line[0] != '#' && line[0] != '@') // comments
                {
                    String contig_name;
                    if (baitlist_format)
                    {
                        vec<String> baitlist_tokens;
                        Tokenize(line, '\t', baitlist_tokens);
                        if (baitlist_tokens.size() < 3)
                        {
                            std::cerr << line << "\n";
                            FatalErr(intervals + " does not appear to have "
                                "proper bait-list formatting");
                        }
                        contig_name = baitlist_tokens[0];
                        init_start.push_back(baitlist_tokens[1].Int() - 1);
                        init_stop.push_back(baitlist_tokens[2].Int());
                    }
                    else
                    {
                        if (line.Contains("\t"))
                        {
                            line = line.Before("\t");
                        }
                        line.GlobalReplaceBy(",", "");
                        contig_name = line.Before(":");
                        init_start.push_back(line.Between(":", "-").Int());
                        init_stop.push_back(line.RevAfter("-").Int());                    
                    }
                    if (init_start.back() > init_stop.back())
                    {
                        std::cerr << line << " is an illegal interval\n";
                        CRD::exit(EXIT_FAILURE);
                    }
                    int pos = Position(names, contig_name);
                    if (pos < 0)
                    {
                        std::cerr << "Failure to process intervals "
                             << intervals << "\n";
                        std::cerr << "unable to locate chromosome "
                            << line.Before(":")
                            << " in reference\n";
                        CRD::exit(EXIT_FAILURE);
                    }
                    
                    if (init_start.back() >= bases[pos].size() ||
                        init_stop.back() > bases[pos].size())
                    {
                        std::cerr << "Interval " << line << " exceeds the "
                            << "boundaries of " << names[pos] << " (size = "
                            << bases[pos].size() << ")\n";
                        CRD::exit(EXIT_FAILURE);
                    }
        
                    init_chr.push_back(pos);
                }
            }
        }
        
        // sort by chromosome
        SortSync(init_chr, init_start, init_stop);
    
        // sort within chromosomes
        for (size_t c = 0; c < init_chr.size(); c++)
        {
            size_t cur_c = init_chr[c];
            vec<size_t> cur_chr;
            vec<size_t> cur_start;
            vec<size_t> cur_stop;
            while (c < init_chr.size() && init_chr[c] == cur_c)
            {
                cur_chr.push_back(init_chr[c]);
                cur_start.push_back(init_start[c]);
                cur_stop.push_back(init_stop[c]);
                c++;
            }
            c--;
            
            SortSync(cur_start, cur_stop);
            
            chr.insert(chr.end(), cur_chr.begin(), cur_chr.end());
            start.insert(start.end(), cur_start.begin(), cur_start.end());
            stop.insert(stop.end(), cur_stop.begin(), cur_stop.end());
        }
        
        return;    
    }
    
    /* Convert a list of intervals to test into the same description as is
       used for typical content motifs. */
    void testintervals2motif(const vec<size_t>& test_chr,
        const vec<size_t>& test_start, const vec<size_t>& test_stop,
        const BC::GenomeReference& gen_ref,
        GenomeInfoCache<bitvector>*& ti_matches,
        GenomeInfoCache<MatchedLocsVec>*& ti_locs)
    {
        const vecbasevector& genome = gen_ref.getBases();
        ti_matches = new GenomeInfoCache<bitvector>(genome, "", true, false,
            true);
        ti_locs = new GenomeInfoCache<MatchedLocsVec>(genome, "", true, false,
            true);
        size_t t = 0;
        for (size_t g = 0; g < genome.size(); g++)
        {
            (*ti_matches)[g].resize(genome[g].size());
            while (t < test_chr.size() && g == test_chr[t])
            {
                (*ti_matches)[g].Set(test_start[t], test_stop[t], true);
                (*ti_locs)[g].push_back(SequenceLoc(test_start[t],
                        test_stop[t] - 1));
                t++;
            }
        }
        return;
    }

    /* Convert a list of intervals to test into a motif describing the
       inverse locations. */
    void testintervals2inversemotif(const vec<size_t>& test_chr,
        const vec<size_t>& test_start, const vec<size_t>& test_stop,
        const BC::GenomeReference& gen_ref,
        GenomeInfoCache<bitvector>*& ti_matches,
        GenomeInfoCache<MatchedLocsVec>*& ti_locs)
    {
        const vecbasevector& genome = gen_ref.getBases();
        ti_matches = new GenomeInfoCache<bitvector>(genome, "", true, false,
            true);
        ti_locs = new GenomeInfoCache<MatchedLocsVec>(genome, "", true, false,
            true);
        size_t t = 0;
        for (size_t g = 0; g < genome.size(); g++)
        {
            (*ti_matches)[g].resize(genome[g].size(), true);
            
            // mask out everything in the intervals
            while (t < test_chr.size() && g == test_chr[t])
            {
                (*ti_matches)[g].Set(test_start[t], test_stop[t], false);
                t++;
            }

            // insert intervals that represent the remainder
            size_t interval_start = 0;
            for (size_t b = 0; b < genome[g].size(); b++)
            {
                if ((*ti_matches)[g][b] && (b + 1 == genome[g].size() || 
                    !(*ti_matches)[g][b + 1]))
                {
                    (*ti_locs)[g].push_back(SequenceLoc(interval_start, b));
                }
                else if (!(*ti_matches)[g][b] && b + 1 < genome[g].size() &&
                    (*ti_matches)[g][b + 1])
                {
                    interval_start = b + 1;
                }
            }
        }
        return;
    }
    
    /* Translate a motif into a filename that is POSIX-safe. */
    String make_cache_name(String motif)
    {
        String cache_name = "";
        for (size_t i = 0; i < motif.size(); i++)
        {
            if (motif[i] == '^')
            {
                cache_name.append("exp");
            }
            else if (motif[i] == '+')
            {
                cache_name.append("plus");
            }
            else if (motif[i] == '=')
            {
                cache_name.append("eq");
            }
            else if (motif[i] == '%')
            {
                cache_name.append("apeq");
            }
            else if (motif[i] == '<')
            {
                cache_name.append("lt");
            }
            else if (motif[i] == '>')
            {
                cache_name.append("gt");
            }
            else if (motif[i] == ')')
            {
                cache_name.append("p_");
            }
            else if (motif[i] == '(')
            {
                cache_name.append("_p");
            }
            else if (motif[i] == '&')
            {
                cache_name.append("_and_");
            }
            else if (motif[i] == '|')
            {
                cache_name.append("_or_");
            }
            else
            {
                cache_name.push_back(motif[i]);
            }
        }
        return cache_name;
    }
    
    /* Initialize the (permanent) motif caches, which contain the locations
       that are matched by these motifs on the specified reference. If caches
       already exist, their paths are returned and no modifications are made. */
    vec<String> initialize_motif_caches(const String& motifmatch_cache_dir,
                                        const String& motifloc_cache_dir,
                                        const vec<String>& motifs,
                                        const BC::GenomeReference& gen_ref,
                                        ostream& out) 
    {
        const vecbasevector& genome = gen_ref.getBases();
        const vecbitvector& refamb = gen_ref.getAmbiguities();
        vec<String> motif_dir_names(motifs.size());
            
        for (size_t mi = 0; mi < motifs.size(); mi++)
        {
            motif_dir_names[mi] = make_cache_name(WhiteSpaceFree(motifs[mi]));
        }
    
        bool motifmatch_cache_dir_exists = IsDirectory(motifmatch_cache_dir);
        bool motifloc_cache_dir_exists = IsDirectory(motifloc_cache_dir);
        if (!motifmatch_cache_dir_exists && !motifloc_cache_dir_exists)
        {
            Mkdir777(motifmatch_cache_dir);
            Mkdir777(motifloc_cache_dir);
        }
        else if (!(motifmatch_cache_dir_exists && motifloc_cache_dir_exists))
        {
            std::cerr << "Either both " << motifmatch_cache_dir << " and "
                 << motifloc_cache_dir << " should exist, or neither should.\n";
            CRD::exit(EXIT_FAILURE);
        }
    
        for (size_t mi = 0; mi < motifs.size(); mi++)
        {
            // this code is just here to fail on unparseable motifs before
            // we make or fetch cache directories for them, defending against
            // picking up caches of invalid partially processed motifs in
            // subsequent runs
            vec<FixedLengthMatch> unused_parse;
            Parse(motifs[mi], unused_parse);
            
            String match_subdir = motifmatch_cache_dir + "/"
                + motif_dir_names[mi];
            String loc_subdir = motifloc_cache_dir + "/"
                + motif_dir_names[mi];
            bool match_exists = IsDirectory(match_subdir);
            bool loc_exists = IsDirectory(loc_subdir);
            
            if (!match_exists && !loc_exists)
            {
                out << "creating cache for " << motifs[mi] << "\n";
                Mkdir777(match_subdir);
                Mkdir777(loc_subdir);
                GenomeInfoCache<bitvector> match_cache(genome, match_subdir,
                                                       false);
                GenomeInfoCache<MatchedLocsVec> loc_cache(genome, loc_subdir,
                                                          false);
                for (size_t chr = 0; chr < genome.size(); chr++)
                {
                    bitvector* chr_match = match_cache.getInfo(chr);
                    MatchedLocsVec* chr_locs = loc_cache.getInfo(chr);
                    chr_match->resize(genome[chr].size(), 0);
                    
                    ComputeMotif(motifs[mi], genome[chr], refamb[chr],
                                 *chr_match, *chr_locs);
                }
                
                std::ofstream finished((match_subdir + 
                    "/finished.txt").c_str());
                finished << motifs[mi] << " cached\n";
                finished.close();
            }
            else if (!(match_exists && loc_exists))
            {
                std::cerr << "Either both " << match_subdir << " and "
                     << loc_subdir << " should exist, or neither should.\n";
                CRD::exit(EXIT_FAILURE);
            }
            else if (!IsRegularFile(match_subdir + "/finished.txt"))
            {
                std::cerr << motifs[mi] << " was not completely initialized in "
                    "the cache, please delete it and try again (or perhaps " 
                    "another instance of BadCoverage is currently creating "
                    "it?)\n";
                CRD::exit(EXIT_FAILURE);
            }
        }
    
        return motif_dir_names;
    }
    
    /* Create a new cache directory, if it does not already exist. */
    String make_cache_dir(String cache_parent_dir, String ref_fasta,
        ostream& out)
    {
        if (!IsSomeSortOfFile(cache_parent_dir))
        {
            out << "creating new cache root directory " << cache_parent_dir
                << "\n";
            Mkdir777(cache_parent_dir.c_str());
        }
        else if (!IsDirectory(cache_parent_dir))
        {
            std::cerr << cache_parent_dir << " is not a directory\n";
            CRD::exit(EXIT_FAILURE);
        }
        
        String new_cache_dir = cache_parent_dir + "/" + ref_fasta;
        
        if (!IsSomeSortOfFile(new_cache_dir))
        {
            out << "creating new cache reference directory " << new_cache_dir
                << "\n";
            Mkpath(new_cache_dir.c_str());
        }
        else if (!IsDirectory(new_cache_dir))
        {
            std::cerr << new_cache_dir << " is not a directory\n";
            CRD::exit(EXIT_FAILURE);
        }
        
        out << "cache = " << new_cache_dir << "\n";
        
        return new_cache_dir;
    }
        
    /* Representation of coverage in a genome location as the output of a
       CoverageIterator is processed. By specifying a bin size > 1, it can
       represent the coverage across non-overlapping continuous regions
       of the genome. */
    class CoverageAccumulator
    {
    public:
        longlong coverage;
        longlong quality;
        longlong mismatches;
        longlong deletions;
        longlong insertions;
        longlong clips;
        longlong skips;
        longlong starts;
        
        CoverageAccumulator(int init_bin_size = 1)
        {
            bin_size = init_bin_size;
            reset();
            return;
        }
        
        CoverageAccumulator& operator+=(const BC::Coverage& cov)
        {
            if (bin_pos == bin_size)
            { 
                reset();
            }
            
            bin_pos++;
            
            coverage += cov.fwd + cov.rev;
            quality += cov.fwd_qual + cov.rev_qual;
            mismatches += cov.fwd_mismatches + cov.rev_mismatches;
            deletions += cov.fwd_deletions + cov.rev_deletions;
            insertions += cov.fwd_insertions + cov.rev_insertions;
            clips += cov.fwd_clips + cov.rev_clips;
            skips += cov.fwd_skips + cov.rev_skips;
            starts += cov.fwd_starts + cov.rev_starts;
            
            if (cov.fwd_qual < 0 || cov.rev_qual < 0)
            {
                quality = -1;
            }
            
            if (cov.fwd_skips < 0 || cov.rev_skips < 0)
            {
                skips = -1;
            }
            
            return *this;
        }
        
    private:
        void reset()
        {
            bin_pos = 0;
            coverage = 0;
            quality = 0;
            mismatches = 0;
            deletions = 0;
            insertions = 0;
            clips = 0;
            skips = 0;
            starts = 0;
            return;
        }
    
        int bin_size;
        int bin_pos;
    };
    
    /* A CoverageMotif determines if a genome location is "in" or "out" based
       on some property of the local coverage. It can have per-base or coarser
       resolution depending on the CoverageAccumulator passed to the inMotif
       method. Motifs based on coverage only become statistically valid
       given high coverage. This super-class manages accumulation of the
       statistics and the intervals of the genome that fit the motif
       definition. */
    class CoverageMotif
    {
    public:
        CoverageMotif(String init_motif_name,
            const vec<String>& init_genome_motifs,
            const vecbitvector& init_refamb,
            const vecbitvector& init_exclusions) :
            motif_name(init_motif_name),
            genome_motifs(init_genome_motifs),
            covmotifsums(0),
            covmotif_nongenomemotif_sums(0),
            covmotif_ingenomemotif_sums(init_genome_motifs.size()),
            genome_N(0),
            in_covmotif(false),
            in_covmotif_nongenomemotif(false),
            refamb(init_refamb),
            exclusions(init_exclusions)
        {
            return;
        }
    
        virtual ~CoverageMotif() {}
    
        void process_coverage(size_t contig, size_t base,
            const CoverageAccumulator& cov,
            const vec<bool>& genome_motif_matched)
        {
            if (!refamb[contig][base] && !exclusions[contig][base])
            {
                if (inMotif(cov))
                {
                    covmotifsums++;
                    
                    if (!in_covmotif || covmotif_intervals.back()[0] != contig)
                    {
                        vec<size_t> new_interval(3, 0);
                        new_interval[0] = contig;
                        new_interval[1] = base;
                        new_interval[2] = base;
                        covmotif_intervals.push_back(new_interval);
                        in_covmotif = true;
                    }
                    ForceAssert(covmotif_intervals.back()[2] == base);
                    covmotif_intervals.back()[2]++;
                    
                    bool any_genome_motif_matched = false;
                    for (size_t mi = 0; mi < genome_motif_matched.size(); mi++)
                    {
                        if (genome_motif_matched[mi])
                        {
                            covmotif_ingenomemotif_sums[mi]++;
                            any_genome_motif_matched = true;
                        }
                    }
                    
                    if (!any_genome_motif_matched)
                    {
                        covmotif_nongenomemotif_sums++;
                        
                        if (!in_covmotif_nongenomemotif || 
                            covmotif_nongenomemotif_intervals.back()[0] != 
                            contig)
                        {
                            vec<size_t> new_interval(3, 0);
                            new_interval[0] = contig;
                            new_interval[1] = base;
                            new_interval[2] = base;
                            covmotif_nongenomemotif_intervals.push_back
                                (new_interval);
                            in_covmotif_nongenomemotif = true;
                        }
                        ForceAssert(covmotif_nongenomemotif_intervals.back()[2] 
                            == base);
                        covmotif_nongenomemotif_intervals.back()[2]++;
                    }
                    else
                    {
                        in_covmotif_nongenomemotif = false;
                    }
                }
                else
                {
                    in_covmotif = false;
                    in_covmotif_nongenomemotif = false;
                }
                
                genome_N++;
            }
            else
            {
                in_covmotif = false;
                in_covmotif_nongenomemotif = false;
            }
            
            return;
        }
        
        virtual bool inMotif(const CoverageAccumulator& cov) = 0;
    
        String name()
        {
            return motif_name;
        }
        
        void print(ostream& out)
        {
            out << std::string(80, '-') << "\n";
            double cm_percent = static_cast<double>(covmotifsums) / 
                genome_N * 100;
            out << motif_name << " fraction  = " << cm_percent << "%\n";
                    
            double cm_ngm_percent =
                static_cast<double>(covmotif_nongenomemotif_sums) / 
                genome_N * 100;
            out << motif_name << " nonmotif fraction  = " << cm_ngm_percent
                << "%\n";
            
            double nongmotif_percent =
                static_cast<double>(covmotif_nongenomemotif_sums) /
                    covmotifsums * 100;
            out << motif_name << " uncovered by any motif = "
                << nongmotif_percent << "%\n";

            for (size_t mi = 0; mi < genome_motifs.size(); mi++)
            {
                double gmotif_percent =
                    static_cast<double>(covmotif_ingenomemotif_sums[mi])
                    / covmotifsums * 100;
                out << motif_name << " covered by " << genome_motifs[mi]
                    << " = " << gmotif_percent << "%\n";
            }

            vec<size_t> covmotif_lengths(covmotif_intervals.size());
            for (size_t c = 0; c < covmotif_lengths.size(); c++)
            {
                covmotif_lengths[c] = (covmotif_intervals[c][2] - 
                    covmotif_intervals[c][1]);
            }
    
            out << motif_name << " N50 = ";
            if (covmotif_lengths.size() > 0)
            {
                Sort(covmotif_lengths);
                out << N50(covmotif_lengths) << "\n";
            }
            else
            {
                out << "[undefined]\n";
            }

            vec<size_t> covmotif_nongenomemotif_lengths
                (covmotif_nongenomemotif_intervals.size());
            for (size_t c = 0; c < covmotif_nongenomemotif_lengths.size(); c++)
            {
                covmotif_nongenomemotif_lengths[c] =
                    (covmotif_nongenomemotif_intervals[c][2] -
                    covmotif_nongenomemotif_intervals[c][1]);
            }
            
            out << motif_name << " non-motif N50 = ";
            if (covmotif_nongenomemotif_lengths.size() > 0)
            {
                Sort(covmotif_nongenomemotif_lengths);
                out << N50(covmotif_nongenomemotif_lengths) << "\n";
            }
            else
            {
                out << "[undefined]\n";
            }
        }
        
        void print_intervals(ostream& out, const vec<String>& contig_names)
        {
            for (size_t c = 0; c < covmotif_intervals.size(); c++)
            {
                out << contig_names[covmotif_intervals[c][0]] << ":"
                    << covmotif_intervals[c][1] << "-"
                    << covmotif_intervals[c][2] << "\n";
            }
            return;
        }
        
        void print_nongenomemotif_intervals(ostream& out,
            const vec<String>& contig_names)
        {
            for (size_t c = 0; c < covmotif_nongenomemotif_intervals.size(); 
                c++)
            {
                out << contig_names[covmotif_nongenomemotif_intervals[c][0]]
                    << ":"
                    << covmotif_nongenomemotif_intervals[c][1] << "-"
                    << covmotif_nongenomemotif_intervals[c][2] << "\n";
            }
            return;        
        }
        
    private:
        String motif_name;
        const vec<String>& genome_motifs;
        longlong covmotifsums;
        longlong covmotif_nongenomemotif_sums;
        vec<longlong> covmotif_ingenomemotif_sums;
        longlong genome_N;
        vec<vec<size_t> > covmotif_intervals;
        vec<vec<size_t> > covmotif_nongenomemotif_intervals;
        bool in_covmotif;
        bool in_covmotif_nongenomemotif;
        const vecbitvector& refamb;
        const vecbitvector& exclusions;
    };
    
    /* Representing regions where coverage is below a specified threshold. */
    class BadCoverageMotif : public CoverageMotif
    {
    public:
        BadCoverageMotif(double init_threshold,
            const vec<String>& init_genome_motifs,
            const vecbitvector& init_refamb,
            const vecbitvector& init_exclusions) :
            CoverageMotif("bad coverage", init_genome_motifs, init_refamb,
                init_exclusions)
        {
            threshold = init_threshold;
            return;
        }
        
        virtual bool inMotif(const CoverageAccumulator& cov)
        {
            return cov.coverage < threshold;
        }
    
    private:
        double threshold;    
    };
    
    /* Reprsenting regions with no coverage. */
    class UncoveredMotif : public CoverageMotif
    {
    public:
        UncoveredMotif(const vec<String>& init_genome_motifs,
            const vecbitvector& init_refamb,
            const vecbitvector& init_exclusions) : 
            CoverageMotif("uncovered", init_genome_motifs,
                init_refamb, init_exclusions)
        {
            return;
        }
        
        virtual bool inMotif(const CoverageAccumulator& cov)
        {
            return cov.coverage == 0;
        }
    };
    
    /* Representing regions with low quality scores. */
    class BadQualityMotif : public CoverageMotif
    {
    public:
        BadQualityMotif(double init_threshold,
            const vec<String>& init_genome_motifs,
            const vecbitvector& init_refamb,
            const vecbitvector& init_exclusions) :
            CoverageMotif("bad quality", init_genome_motifs,
                init_refamb, init_exclusions)
        {
            threshold = init_threshold;
            return;
        }
        
        virtual bool inMotif(const CoverageAccumulator& cov)
        {
            return cov.quality < threshold * cov.coverage;
        }
    
    private:
        double threshold;
    };
    
    /* Representing regions with a high mismatch rate. */
    class BadMismatchesMotif : public CoverageMotif
    {
    public:
        BadMismatchesMotif(double init_threshold,
            const vec<String>& init_genome_motifs,
            const vecbitvector& init_refamb,
            const vecbitvector& init_exclusions) :
            CoverageMotif("bad mismatches", init_genome_motifs,
                init_refamb, init_exclusions)
        {
            threshold = init_threshold;
            return;
        }
        
        virtual bool inMotif(const CoverageAccumulator& cov)
        {
            return cov.coverage - cov.mismatches < (1 - threshold) * 
                cov.coverage;
        }
    
    private:
        double threshold;
    };
    
    /* Representing regions with a high deletion rate. */
    class BadDeletionsMotif : public CoverageMotif
    {
    public:
        BadDeletionsMotif(double init_threshold,
            const vec<String>& init_genome_motifs,
            const vecbitvector& init_refamb,
            const vecbitvector& init_exclusions) :
            CoverageMotif("bad deletions", init_genome_motifs,
                init_refamb, init_exclusions)
        {
            threshold = init_threshold;
            return;
        }
        
        virtual bool inMotif(const CoverageAccumulator& cov)
        {
            return cov.coverage - cov.deletions < (1 - threshold) * 
                cov.coverage;
        }
    
    private:
        double threshold;
    };
    
    /* Representing regions with a high insertion rate. */
    class BadInsertionsMotif : public CoverageMotif
    {
    public:
        BadInsertionsMotif(double init_threshold,
            const vec<String>& init_genome_motifs,
            const vecbitvector& init_refamb,
            const vecbitvector& init_exclusions) :
            CoverageMotif("bad insertions", init_genome_motifs,
                init_refamb, init_exclusions)
        {
            threshold = init_threshold;
            return;
        }
    
        virtual bool inMotif(const CoverageAccumulator& cov)
        {
            return cov.coverage - cov.insertions < (1 - threshold) * 
                cov.coverage;
        }
    
    private:
        double threshold;
    };
    
    /* Representing regions with a high clipping rate. */
    class BadClipsMotif : public CoverageMotif
    {
    public:
        BadClipsMotif(double init_threshold,
            const vec<String>& init_genome_motifs,
            const vecbitvector& init_refamb,
            const vecbitvector& init_exclusions) :
            CoverageMotif("bad clips", init_genome_motifs, init_refamb,
                init_exclusions)
        {
            threshold = init_threshold;
            return;
        }
        
        virtual bool inMotif(const CoverageAccumulator& cov)
        {
            return cov.coverage - cov.clips < (1 - threshold) * cov.coverage;
        }
    
    private:
        double threshold;
    };
    
    /* Representing regions that are the inverse of another CoverageMotif. */
    class InverseMotif : public CoverageMotif
    {
    public:
        InverseMotif(CoverageMotif& init_source_motif, String init_motif_name,
            const vec<String>& init_genome_motifs,
            const vecbitvector& init_refamb,
            const vecbitvector& init_exclusions) : 
            CoverageMotif(init_motif_name, init_genome_motifs,
                init_refamb, init_exclusions), 
            source_motif(init_source_motif)
        {
            return;
        }
        
        virtual bool inMotif(const CoverageAccumulator& cov)
        {
            return !source_motif.inMotif(cov);
        }
        
    private:
        CoverageMotif& source_motif;
    };
    
    /* Representing regions that are the logical AND of two CoverageMotifs. */
    class IntersectionMotif : public CoverageMotif
    {
    public:
        IntersectionMotif(CoverageMotif& init_source_motifA,
            CoverageMotif& init_source_motifB,
            const vec<String>& init_genome_motifs,
            const vecbitvector& init_refamb,
            const vecbitvector& init_exclusions) :
            CoverageMotif(init_source_motifA.name() + " and "
                + init_source_motifB.name(), init_genome_motifs,
                init_refamb, init_exclusions),
            source_motifA(init_source_motifA),
            source_motifB(init_source_motifB)
        {
            return;        
        }
        
        virtual bool inMotif(const CoverageAccumulator& cov)
        {
            return source_motifA.inMotif(cov) && source_motifB.inMotif(cov);
        }
    
    private:
        CoverageMotif& source_motifA;
        CoverageMotif& source_motifB;
    };
    
    /* Class for accumulating statistics for the GC-bias table. */
    class GCBiasTable
    {
    public:
        GCBiasTable(unsigned int init_winsize, const vecbasevector& init_genome,
            const vecbitvector& init_refamb,
            const vecbitvector& init_exclusions) :
            winsize(init_winsize),
            genome(init_genome),
            refamb(init_refamb),
            exclusions(init_exclusions),
            hundred_gc_windows(101, 0),
            hundred_gc_coverage(101, 0), 
            hundred_gc_mismatches(101, 0),
            hundred_gc_insertions(101, 0), 
            hundred_gc_deletions(101, 0),
            hundred_gc_clips(101, 0), 
            hundred_gc_skips(101, 0),
            hundred_gc_quality(101, 0), 
            hundred_gc_starts(101, 0),
            gc_bases_sum(0),
            gc_window_valid_sum(0),
            gc_window_coverage_sum(0), 
            gc_window_mismatches_sum(0),
            gc_window_insertions_sum(0), 
            gc_window_deletions_sum(0),
            gc_window_clips_sum(0),
            gc_window_skips_sum(0),
            gc_window_quality_sum(0),
            gc_window_starts_sum(0),
            missing_quality(false),
            missing_skips(false)
        {
            return;
        }
        
        void process_coverage(size_t contig, size_t base,
            const CoverageAccumulator& cov)
        {
            if (base == 0)
            {
                gc_bases.clear();
                gc_bases_sum = 0;
                gc_window_valid.clear();
                gc_window_valid_sum = 0;
                gc_window_coverage.clear();
                gc_window_coverage_sum = 0;
                gc_window_mismatches.clear();
                gc_window_mismatches_sum = 0;
                gc_window_insertions.clear();
                gc_window_insertions_sum = 0;
                gc_window_deletions.clear();
                gc_window_deletions_sum = 0;
                gc_window_clips.clear();
                gc_window_clips_sum = 0;
                gc_window_skips.clear();
                gc_window_skips_sum = 0;
                gc_window_quality.clear();
                gc_window_quality_sum = 0;
                gc_window_starts.clear();
                gc_window_starts_sum = 0;
            }
            
            if (!refamb[contig][base] && !exclusions[contig][base])
            {
                gc_window_valid.push_back(true);
                gc_window_valid_sum++;
                gc_bases.push_back(IsGC(genome[contig][base]));
                gc_bases_sum += gc_bases.back();
                gc_window_coverage.push_back(cov.coverage);
                gc_window_coverage_sum += cov.coverage;
                gc_window_mismatches.push_back(cov.mismatches);
                gc_window_mismatches_sum += cov.mismatches;
                gc_window_insertions.push_back(cov.insertions);
                gc_window_insertions_sum += cov.insertions;
                gc_window_deletions.push_back(cov.deletions);
                gc_window_deletions_sum += cov.deletions;
                gc_window_clips.push_back(cov.clips);
                gc_window_clips_sum += cov.clips;
                gc_window_skips.push_back(cov.skips);
                gc_window_skips_sum += cov.skips;
                gc_window_quality.push_back(cov.quality);
                gc_window_quality_sum += cov.quality;
                gc_window_starts.push_back(cov.starts);
                gc_window_starts_sum += cov.starts;
                
                if (cov.quality < 0)
                {
                    missing_quality = true;
                }
                
                if (cov.skips < 0)
                {
                    missing_skips = true;
                }
            }
            else
            {
                gc_window_valid.push_back(false);
                gc_bases.push_back(false);
                gc_window_coverage.push_back(0);
                gc_window_mismatches.push_back(0);
                gc_window_insertions.push_back(0);
                gc_window_deletions.push_back(0);
                gc_window_clips.push_back(0);
                gc_window_skips.push_back(0);
                gc_window_quality.push_back(0);
                gc_window_starts.push_back(0);
            }
            
            if (gc_window_valid.size() > winsize)
            {
                gc_window_valid_sum -= gc_window_valid.front();
                gc_window_valid.pop_front();
                gc_bases_sum -= gc_bases.front();
                gc_bases.pop_front();
                gc_window_coverage_sum -= gc_window_coverage.front();
                gc_window_coverage.pop_front();
                gc_window_mismatches_sum -= gc_window_mismatches.front();
                gc_window_mismatches.pop_front();
                gc_window_insertions_sum -= gc_window_insertions.front();
                gc_window_insertions.pop_front();
                gc_window_deletions_sum -= gc_window_deletions.front();
                gc_window_deletions.pop_front();
                gc_window_clips_sum -= gc_window_clips.front();
                gc_window_clips.pop_front();
                gc_window_skips_sum -= gc_window_skips.front();
                gc_window_skips.pop_front();
                gc_window_quality_sum -= gc_window_quality.front();
                gc_window_quality.pop_front();
                gc_window_starts_sum -= gc_window_starts.front();
                gc_window_starts.pop_front();                    
            }
            
            ForceAssert(gc_window_valid.size() <= winsize);
            
            if (winsize > 0 && gc_window_valid_sum == winsize)
            {
                int gc_index = 100 * gc_bases_sum / winsize;
                hundred_gc_windows[gc_index]++;
                hundred_gc_coverage[gc_index] += gc_window_coverage_sum;
                hundred_gc_mismatches[gc_index] += gc_window_mismatches_sum;
                hundred_gc_insertions[gc_index] += gc_window_insertions_sum;
                hundred_gc_deletions[gc_index] += gc_window_deletions_sum;
                hundred_gc_clips[gc_index] += gc_window_clips_sum;
                hundred_gc_skips[gc_index] += gc_window_skips_sum;
                hundred_gc_quality[gc_index] += gc_window_quality_sum;
                hundred_gc_starts[gc_index] += gc_window_starts_sum;
            }            
            return;
        }
        
        void print(ostream& out)
        {
            double mean_rel_coverage = 0;
            double mean_rel_mismatches = 0;
            double mean_rel_insertions = 0;
            double mean_rel_deletions = 0;
            double mean_rel_clips = 0;
            double mean_rel_skips = 0;
            double mean_rel_quality = 0;
            double mean_rel_starts = 0;
            longlong total_windows = 0;
            longlong total_coverage = 0;
            for (size_t gc = 0; gc < hundred_gc_coverage.size(); gc++)
            {
                mean_rel_coverage += 
                    static_cast<double>(hundred_gc_coverage[gc]);
                mean_rel_mismatches +=
                    static_cast<double>(hundred_gc_mismatches[gc]);
                mean_rel_insertions +=
                    static_cast<double>(hundred_gc_insertions[gc]);
                mean_rel_deletions += 
                    static_cast<double>(hundred_gc_deletions[gc]);
                mean_rel_clips += static_cast<double>(hundred_gc_clips[gc]);
                mean_rel_skips += static_cast<double>(hundred_gc_skips[gc]);
                mean_rel_quality += static_cast<double>(hundred_gc_quality[gc]);
                mean_rel_starts += static_cast<double>(hundred_gc_starts[gc]);
                total_windows += hundred_gc_windows[gc];
                total_coverage += hundred_gc_coverage[gc];
            }
            mean_rel_coverage /= total_windows;
            mean_rel_mismatches /= total_coverage;
            mean_rel_insertions /= total_coverage;
            mean_rel_deletions /= total_coverage;
            mean_rel_clips /= total_coverage;
            mean_rel_skips /= total_coverage;
            mean_rel_quality /= total_coverage;
            mean_rel_starts /= total_windows;
            
            if (missing_quality)
            {
                mean_rel_quality = numeric_limits<double>::quiet_NaN();
            }
            
            if (missing_skips)
            {
                mean_rel_skips = numeric_limits<double>::quiet_NaN();
            }
            
            out << std::string(80, '-') << "\n";
            out << "GC%\tWIN\tCOV\tMIS\tINS\tDEL\tCLP\tSKP\tQUA\tRST\n";
            for (size_t gc = 0; gc <= 100; gc++)
            {
                if (hundred_gc_windows[gc] > 0)
                {
                    double rel_windows = 
                        static_cast<double>(hundred_gc_windows[gc])
                        / total_windows;
                    double rel_coverage =
                        static_cast<double>(hundred_gc_coverage[gc])
                        / hundred_gc_windows[gc] / mean_rel_coverage;
                    double rel_mismatches =
                        static_cast<double>(hundred_gc_mismatches[gc])
                        / hundred_gc_coverage[gc] / mean_rel_mismatches;
                    double rel_insertions =
                        static_cast<double>(hundred_gc_insertions[gc])
                        / hundred_gc_coverage[gc] / mean_rel_insertions;
                    double rel_deletions =
                        static_cast<double>(hundred_gc_deletions[gc])
                        / hundred_gc_coverage[gc] / mean_rel_deletions;
                    double rel_clips = static_cast<double>(hundred_gc_clips[gc])
                        / hundred_gc_coverage[gc] / mean_rel_clips;
                    double rel_skips = static_cast<double>(hundred_gc_skips[gc])
                        / hundred_gc_coverage[gc] / mean_rel_skips;
                    double rel_quality = 
                        static_cast<double>(hundred_gc_quality[gc])
                        / hundred_gc_coverage[gc] / mean_rel_quality;
                    double rel_starts = 
                        static_cast<double>(hundred_gc_starts[gc])
                        / hundred_gc_windows[gc] / mean_rel_starts;
                    out << gc << "\t"
                        << rel_windows << "\t"
                        << rel_coverage << "\t"
                        << rel_mismatches << "\t"
                        << rel_insertions << "\t"
                        << rel_deletions << "\t"
                        << rel_clips << "\t"
                        << rel_skips << "\t"
                        << rel_quality << "\t"
                        << rel_starts << "\n";
                }
            }
            out << "TOTAL " << winsize << "-BASE WINDOWS = "
                << total_windows << "\n";
            out << "TOTAL " << winsize << "-BASE COVERAGE = "
                << total_coverage << "\n";
            out << "MEAN " << winsize << "-BASE COVERAGE = "
                << mean_rel_coverage << "\n";
            out << "MEAN " << winsize << "-BASE MISMATCHES = " 
                << mean_rel_mismatches << "\n";
            out << "MEAN " << winsize << "-BASE INSERTIONS = "
                << mean_rel_insertions << "\n";
            out << "MEAN " << winsize << "-BASE DELETIONS = "
                << mean_rel_deletions << "\n";
            out << "MEAN " << winsize << "-BASE CLIPS = "
                << mean_rel_clips << "\n";
            out << "MEAN " << winsize << "-BASE SKIPS = "
                << mean_rel_skips << "\n";
            out << "MEAN " << winsize << "-BASE QUALITY = "
                << mean_rel_quality << "\n";
            out << "MEAN " << winsize << "-BASE READ STARTS = "
                << mean_rel_starts << "\n";
        }
        
    private:
        unsigned int winsize;
        const vecbasevector& genome;
        const vecbitvector& refamb;
        const vecbitvector& exclusions;
        vec<longlong> hundred_gc_windows;
        vec<longlong> hundred_gc_coverage;
        vec<longlong> hundred_gc_mismatches;
        vec<longlong> hundred_gc_insertions;
        vec<longlong> hundred_gc_deletions;
        vec<longlong> hundred_gc_clips;
        vec<longlong> hundred_gc_skips;
        vec<longlong> hundred_gc_quality;
        vec<longlong> hundred_gc_starts;
        std::deque<bool> gc_bases;
        int gc_bases_sum;
        std::deque<bool> gc_window_valid;
        unsigned int gc_window_valid_sum;
        std::deque<int> gc_window_coverage;
        int gc_window_coverage_sum;
        std::deque<int> gc_window_mismatches;
        int gc_window_mismatches_sum;
        std::deque<int> gc_window_insertions;
        int gc_window_insertions_sum;
        std::deque<int> gc_window_deletions;
        int gc_window_deletions_sum;
        std::deque<int> gc_window_clips;
        int gc_window_clips_sum;
        std::deque<int> gc_window_skips;
        int gc_window_skips_sum;
        std::deque<int> gc_window_quality;
        int gc_window_quality_sum;
        std::deque<int> gc_window_starts;
        int gc_window_starts_sum;
        bool missing_quality;
        bool missing_skips;
    };
    
    /* Class for accumulating statistics for the read-start latent variance
       analysis. */
    class LatentVarianceAnalysis
    {
    public:
        LatentVarianceAnalysis(bool init_report_latent_variance,
            bool init_report_latent_variance_ci,
            const vecbitvector& init_refamb,
            const vecbitvector& init_exclusions) :
            report_latent_variance(init_report_latent_variance),
            report_latent_variance_ci(init_report_latent_variance_ci),
            read_starts_mean(0),
            read_starts_variance(0),
            read_starts_M2(0),
            total_read_starts(0),
            read_starts_kurtosis(0),
            read_starts_kurtosis_M2(0),
            read_starts_kurtosis_M3(0),
            read_starts_kurtosis_M4(0),
            N(0),
            refamb(init_refamb),
            exclusions(init_exclusions)
        {
            return;
        }
        
        void process_coverage(size_t contig, size_t base,
            const CoverageAccumulator& cov)
        {
            // on-line mean/variance calculation, originally from
            // Knuth/Welford (see http://en.wikipedia.org/wiki/
            // Algorithms_for_calculating_variance)

            // on-line kurtosis added from (http://en.wikipedia.org/
            // wiki/Algorithms_for_calculating_variance
            // #Higher-order_statistics)
            // note that N hasn't been incremented yet
            if (!refamb[contig][base] && !exclusions[contig][base])
            {
                double delta = cov.starts - read_starts_mean;
                double delta_N = delta / (N + 1);
                double delta_N2 = delta_N * delta_N;
                double term1 = delta * delta_N * N;
                read_starts_mean += delta_N;
                read_starts_M2 += delta * (cov.starts - read_starts_mean);
                read_starts_variance = read_starts_M2 / N;
                total_read_starts += cov.starts;
    
                read_starts_kurtosis_M4 +=
                    term1 * delta_N2 * ((N + 1) * (N + 1) - 3 * (N + 1) + 3)
                    + 6 * delta_N2 * read_starts_kurtosis_M2
                    - 4 * delta_N * read_starts_kurtosis_M3;
                read_starts_kurtosis_M3 +=
                    term1 * delta_N * (N + 1 - 2)
                    - 3 * delta_N * read_starts_kurtosis_M2;
                read_starts_kurtosis_M2 += term1;
                read_starts_kurtosis = ((N + 1) * read_starts_kurtosis_M4) /
                    (read_starts_kurtosis_M2 * read_starts_kurtosis_M2) - 3;
                    
                N++;
            }
            
            return;
        }
        
        void print(ostream& out)
        {
            out << "var read starts = " << read_starts_variance << "\n";
            if (report_latent_variance)
            {
                double latent_placement_var =
                    ((((read_starts_variance / pow(read_starts_mean, 2))
                       - (1 / read_starts_mean))
                      + 1 / (static_cast<double>(total_read_starts)))
                     * (static_cast<double>(total_read_starts)
                        / static_cast<double>(total_read_starts - 1)));
                
                double rs_var_var = pow(read_starts_variance, 2)
                    * (2 / (N - 1) + read_starts_kurtosis / N);
            
                double latent_placement_var_var =
                    pow(total_read_starts, 2) / pow(total_read_starts - 1, 2)
                    * pow(read_starts_mean, -4) * rs_var_var;
            
                double latent_placement_var_ci = 1.96 * 
                    sqrt(latent_placement_var_var);
        
                out << "latent read-placement probability variance = "
                    << latent_placement_var;
                if (report_latent_variance_ci)
                {
                    out << " [+/- " << latent_placement_var_ci << "]";
                }
                out << "\n";
            }
            return;
        }
        
    private:
        bool report_latent_variance;
        bool report_latent_variance_ci;
        double read_starts_mean;
        double read_starts_variance;
        double read_starts_M2;
        longlong total_read_starts;
    
        double read_starts_kurtosis;
        double read_starts_kurtosis_M2;
        double read_starts_kurtosis_M3;
        double read_starts_kurtosis_M4;
        longlong N;
        
        const vecbitvector& refamb;
        const vecbitvector& exclusions;
    };
    
    /* Class for accumulating statistics across the entire genome. */
    class BulkStatistics
    {
    public:
        BulkStatistics(const vecbitvector& init_refamb,
            const vecbitvector& init_exclusions) :
            refamb(init_refamb),
            exclusions(init_exclusions),
            total(0),
            total_sq(0),
            total_quality(0),
            total_mismatches(0),
            total_deletions(0),
            total_insertions(0),
            total_clips(0),
            total_skips(0),
            total_starts(0),
            missing_quality_flag(false),
            missing_skips_flag(false),
            num_unamb_bases(0),
            num_ref_bases(0),
            num_excluded_bases(0)

        {
            return;
        }
        
        void process_coverage(size_t contig, size_t base,
            const CoverageAccumulator& cov)
        {
            num_ref_bases++;
            if (exclusions[contig][base] && !refamb[contig][base])
            {
                num_excluded_bases++;
            }
            
            if (!refamb[contig][base])
            {
                num_unamb_bases++;
                if (!exclusions[contig][base])
                {
                    if (cov.quality < 0)
                    {
                        missing_quality_flag = true;
                    }
                    
                    if (cov.skips < 0)
                    {
                        missing_skips_flag = true;
                    }
                    
                    total += cov.coverage;
                    total_sq += cov.coverage * cov.coverage;
                    total_quality += cov.quality;
                    total_mismatches += cov.mismatches;
                    total_deletions += cov.deletions;
                    total_insertions += cov.insertions;
                    total_clips += cov.clips;
                    total_skips += cov.skips;
                    total_starts += cov.starts;
                }
            }
            return;
        }
        
        double mean_coverage()
        {
            return double(total) / double(num_unamb_bases - num_excluded_bases);
        }

        double mean_quality()
        {
            return double(total_quality) / double(total);
        }

        double mean_mismatches()
        {
            return double(total_mismatches) / double(total);
        }
        
        double mean_deletions()
        {
            return double(total_deletions) / double(total);
        }
        
        double mean_insertions()
        {
            return double(total_insertions) / double(total);
        }
        
        double mean_clips()
        {
            return double(total_clips) / double(total);
        }
        
        double mean_skips()
        {
            return double(total_skips) / double(total);
        }
        
        double mean_starts()
        {
            return double(total_starts) /
                double(num_unamb_bases - num_excluded_bases);
        }
        
        bool missing_quality()
        {
            return missing_quality_flag;
        }
        
        bool missing_skips()
        {
            return missing_skips_flag;
        }
        
        void print(ostream& out)
        { 
            double var_cov = total_sq / double(num_unamb_bases - 
                num_excluded_bases) - pow(mean_coverage(), 2);
            out << "number of unambiguous reference bases = "
                << num_unamb_bases << "\n";
            out << "number of excluded reference bases = "
                << num_excluded_bases << "\n";
            out << "number of included reference bases = "
                << (num_unamb_bases - num_excluded_bases) << "\n";
            out << "total mapped bases = " << total << "\n";
            longlong total_read_bases = total + total_insertions + total_clips;
            out << "coverage per read base = "
                << total / double(total_read_bases) << "\n";
            out << "mismatches per read base = "
                << total_mismatches / double(total_read_bases) << "\n";
            out << "deletions per read base = "
                << total_deletions / double(total_read_bases) << "\n";
            out << "insertions per read base = "
                << total_insertions / double(total_read_bases) << "\n";
            out << "clips per read base = "
                << total_clips / double(total_read_bases) << "\n";
            out << "skips per read base = ";
            if (missing_skips())
            {
                out << "N/A\n";
            }
            else
            {
                out << total_skips / double(total_read_bases) << "\n";
            }
            out << "mean coverage = " << mean_coverage() << "\n";
            out << "var coverage = " << var_cov << "\n";
            out << "mean quality = ";
            if (missing_quality())
            {
                out << "N/A";
            }
            else
            {
                out << mean_quality();
            }
            out << "\n";
            out << "mean mismatches per coverage = " << mean_mismatches()
                << "\n";
            out << "mean deletions per coverage = " << mean_deletions() << "\n";
            out << "mean insertions per coverage = " << mean_insertions()
                << "\n";
            out << "mean clips per coverage = " << mean_clips() << "\n";
            out << "mean skips per coverage = ";
            if (missing_skips())
            {
                out << "N/A\n";
            }
            else
            {
                out << mean_skips() << "\n";
            }
            out << "mean read starts = " << mean_starts() << "\n";

            return;
        }
        
    private:
        const vecbitvector& refamb;
        const vecbitvector& exclusions;
        longlong total;
        longlong total_sq;
        longlong total_quality;
        longlong total_mismatches;
        longlong total_deletions;
        longlong total_insertions;
        longlong total_clips;
        longlong total_skips;
        longlong total_starts;
        bool missing_quality_flag;
        bool missing_skips_flag;
        longlong num_unamb_bases;
        longlong num_ref_bases;
        longlong num_excluded_bases;
    };
    
    /* A class for accumulating the statistics of a particular genome context,
       e.g. "bases in 100-base regions that have less than 10% G+C". */
    class ContentMotifStatistics
    {
    public:
        ContentMotifStatistics(String init_motif, size_t init_motif_number,
            GenomeInfoCache<bitvector>& init_match_cache, 
            GenomeInfoCache<MatchedLocsVec>& init_loc_cache,
            const vecbitvector& init_refamb,
            const vecbitvector& init_exclusions,
            const vec<String>& init_contig_names,
            const bool init_split_intervals) :
            motif(init_motif),
            motif_number(init_motif_number),
            match_cache(init_match_cache),
            loc_cache(init_loc_cache),
            refamb(init_refamb),
            exclusions(init_exclusions),
            contig_names(init_contig_names),
            split_intervals(init_split_intervals),
            current_contig(numeric_limits<size_t>::max()),
            motif_cov_total(0),
            motif_qual_total(0),
            motif_mismatch_total(0),
            motif_deletion_total(0),
            motif_insertion_total(0),
            motif_clip_total(0),
            motif_skip_total(0),
            motif_start_total(0),
            motif_N(0),
            interval_cov_total(0),
            interval_qual_total(0),
            interval_mismatch_total(0),
            interval_deletion_total(0),
            interval_insertion_total(0),
            interval_clip_total(0),
            interval_skip_total(0),
            interval_start_total(0),
            interval_N(0),
            genome_cov_total(0),
            genome_qual_total(0),
            genome_mismatch_total(0),
            genome_deletion_total(0),
            genome_insertion_total(0),
            genome_clip_total(0),
            genome_skip_total(0),
            genome_start_total(0),
            genome_N(0),
            total_motif_matches(0),
            motif_var_size(false),
            ms_loc(0),
            missing_quality(false),
            missing_skips(false)
        {
            return;
        }
        
        bool process_coverage(size_t contig, size_t base,
            const CoverageAccumulator& cov)
        {
            if (contig != current_contig)
            {
                in_motif = match_cache.getInfo(contig);
                matched_locs = loc_cache.getInfo(contig);
                valid_locs.resize(matched_locs->size());
                valid_locs.Set(0, valid_locs.size(), 1);

                if (contig >= interval_cov_total.size())
                {
                    interval_cov_total.resize(contig + 1);
                    interval_qual_total.resize(contig + 1);
                    interval_mismatch_total.resize(contig + 1);
                    interval_deletion_total.resize(contig + 1);
                    interval_insertion_total.resize(contig + 1);
                    interval_clip_total.resize(contig + 1);
                    interval_skip_total.resize(contig + 1);
                    interval_start_total.resize(contig + 1);
                    interval_N.resize(contig + 1);
                }
                
                if (split_intervals && interval_cov_total[contig].size() == 0)
                {
                    interval_cov_total[contig].resize(matched_locs->size(), 0);
                    interval_qual_total[contig].resize(matched_locs->size(), 0);
                    interval_mismatch_total[contig].resize(matched_locs->size(),
                        0);
                    interval_deletion_total[contig].resize(matched_locs->size(), 
                        0);
                    interval_insertion_total[contig].resize
                        (matched_locs->size(), 0);
                    interval_clip_total[contig].resize(matched_locs->size(), 0);
                    interval_skip_total[contig].resize(matched_locs->size(), 0);
                    interval_start_total[contig].resize(matched_locs->size(), 
                        0);
                    interval_N[contig].resize(matched_locs->size(), 0);
                }
                
                for (size_t k = 0; k < matched_locs->size(); k++)
                {
                    size_t match_size = (*matched_locs)[k].second
                        - (*matched_locs)[k].first + 1;
                    if (motif_coverage_counts.size() == 0)
                    {
                        motif_coverage_counts.resize(match_size, 0);
                        motif_quality_counts.resize(match_size, 0);
                    }
                    
                    if (!motif_var_size &&
                        motif_coverage_counts.size() != match_size)
                    {
                        motif_var_size = true;
                    }
                    
                    if (k > 0 && (*matched_locs)[k].first <
                        (*matched_locs)[k - 1].first)
                    {
                        std::cerr << "motif match locations are unsorted\n";
                        CRD::exit(EXIT_FAILURE);
                    }
                    
                    // mark motif locations that cross ambiguous or excluded
                    // regions
                    for (size_t i = (*matched_locs)[k].first;
                        i <= (*matched_locs)[k].second && valid_locs[k];
                        i++)
                    {
                        if (refamb[contig][i] || exclusions[contig][i])
                        {
                            valid_locs.Set(k, 0);
                        }
                    }
                                        
                    if (valid_locs[k])
                    {
                        total_motif_matches++;
                    }
                }
    
                ms_loc = 0;

                current_contig = contig;
            }
            
            bool genome_motif_matched = false;

            if (!refamb[contig][base] && !exclusions[contig][base])
            {
                if (cov.quality < 0)
                {
                    missing_quality = true;
                }
                
                if (cov.skips < 0)
                {
                    missing_skips = true;
                } 
                
                genome_cov_total += cov.coverage;
                genome_qual_total += cov.quality;
                genome_mismatch_total += cov.mismatches;
                genome_deletion_total += cov.deletions;
                genome_insertion_total += cov.insertions;
                genome_clip_total += cov.clips;
                genome_skip_total += cov.skips;
                genome_start_total += cov.starts;
                genome_N++;
    
                if ((*in_motif)[base])
                {
                    genome_motif_matched = true;
                    motif_cov_total += cov.coverage;
                    motif_qual_total += cov.quality;
                    motif_mismatch_total += cov.mismatches;
                    motif_deletion_total += cov.deletions;
                    motif_insertion_total += cov.insertions;
                    motif_clip_total += cov.clips;
                    motif_skip_total += cov.skips;
                    motif_start_total += cov.starts;
                    motif_N++;
                }
    
                // discard motif locations that we've already passed
                while (ms_loc < matched_locs->size() &&
                    static_cast<longlong>(base) >
                    (*matched_locs)[ms_loc].first &&
                    static_cast<longlong>(base) > 
                    (*matched_locs)[ms_loc].second)
                {
                    ms_loc++;
                }
                
                // add coverage from any relevant locations
                for (size_t k = ms_loc; k < matched_locs->size() && 
                    static_cast<longlong>(base) >= (*matched_locs)[k].first;
                    k++)
                {
                    if (base >= (*matched_locs)[k].first &&
                        base <= (*matched_locs)[k].second && valid_locs[k])
                    {
                        if (!motif_var_size)
                        {
                            int match_ind = base - (*matched_locs)[k].first;
                            motif_coverage_counts[match_ind] += cov.coverage;
                            motif_quality_counts[match_ind] += cov.quality;
                        }
                        if (split_intervals)
                        {
                            interval_cov_total[contig][k] += cov.coverage;
                            interval_qual_total[contig][k] += cov.quality;
                            interval_mismatch_total[contig][k] += 
                                cov.mismatches;
                            interval_deletion_total[contig][k] += cov.deletions;
                            interval_insertion_total[contig][k] += 
                                cov.insertions;
                            interval_clip_total[contig][k] += cov.clips;
                            interval_skip_total[contig][k] += cov.skips;
                            interval_start_total[contig][k] += cov.starts;
                            interval_N[contig][k]++;
                        }
                    }
                }
            }
                       
            return genome_motif_matched;
        }
        
        static void print_table_header(ostream& out)
        {
            out << "motif number\t" << "motif\t" << "genome coverage\t"
                << "bases\t" << "locations\t" << "relative coverage\t"
                << "relative quality\t" << "relative mismatches\t"
                << "relative deletions\t" << "relative insertions\t"
                << "relative clips\t" << "relative skips\t"
                << "relative starts\n";
        }
        
        void print_table_data(ostream& out)
        {
            double mean_cov = static_cast<double>(genome_cov_total) /
                genome_N;
            double mean_qual = static_cast<double>(genome_qual_total) / 
                genome_cov_total;
            double mean_mismatch = static_cast<double>(genome_mismatch_total) / 
                genome_cov_total;
            double mean_deletions = static_cast<double>(genome_deletion_total) / 
                genome_cov_total;
            double mean_insertions = static_cast<double>(genome_insertion_total) 
                / genome_cov_total;
            double mean_clips = static_cast<double>(genome_clip_total) / 
                genome_cov_total;
            double mean_skips = static_cast<double>(genome_skip_total) /
                genome_cov_total;
            double mean_starts = static_cast<double>(genome_start_total) /
                genome_N;
   
            double motif_rel_coverage = static_cast<double>(motif_cov_total) / 
                motif_N / mean_cov;
            double motif_rel_quality = static_cast<double>(motif_qual_total) /
                motif_cov_total / mean_qual;
            double motif_rel_mis = static_cast<double>(motif_mismatch_total) /
                motif_cov_total / mean_mismatch;
            double motif_rel_del = static_cast<double>(motif_deletion_total) /
                motif_cov_total / mean_deletions;
            double motif_rel_ins = static_cast<double>(motif_insertion_total) /
                motif_cov_total / mean_insertions;
            double motif_rel_clp = static_cast<double>(motif_clip_total) /
                motif_cov_total / mean_clips;
            double motif_rel_skp = static_cast<double>(motif_skip_total) /
                motif_cov_total / mean_skips;
            double motif_rel_start = static_cast<double>(motif_start_total) /
                motif_N / mean_starts;
            
            out << motif_number << "\t"
                << motif << "\t"
                << (static_cast<double>(motif_N) / genome_N * 100) << "%\t"
                << motif_N << "\t"
                << total_motif_matches << "\t"
                << motif_rel_coverage << "\t";
            if (missing_quality)
            {
                out << "N/A\t";
            }
            else
            {
                out << motif_rel_quality << "\t";
            }
            out << motif_rel_mis << "\t"
                << motif_rel_del << "\t"
                << motif_rel_ins << "\t"
                << motif_rel_clp << "\t";
            if (missing_skips)
            {
                out << "N/A\t";
            }
            else
            {
                out << motif_rel_skp << "\t";
            }
            out << motif_rel_start << "\n";

            for (size_t c = 0; split_intervals && c < interval_cov_total.size(); 
                c++)
            {
                matched_locs = loc_cache.getInfo(c);
                for (size_t i = 0; i < interval_cov_total[c].size(); i++)
                {
                    size_t interval_start = (*matched_locs)[i].first;
                    size_t interval_end = (*matched_locs)[i].second;
                    bool interval_valid = true;
                    
                    for (size_t b = interval_start; b <= interval_end; b++)
                    {
                        if (refamb[c][b] || exclusions[c][b])
                        {
                            interval_valid = false;
                        }
                    }

                    if (interval_valid && interval_N[c][i] > 0)
                    {
                        double interval_rel_coverage = 
                            static_cast<double>(interval_cov_total[c][i]) / 
                            interval_N[c][i] / mean_cov;
                        double interval_rel_quality = 
                            static_cast<double>(interval_qual_total[c][i]) /
                            interval_cov_total[c][i] / mean_qual;
                        double interval_rel_mis = 
                            static_cast<double>(interval_mismatch_total[c][i]) /
                            interval_cov_total[c][i] / mean_mismatch;
                        double interval_rel_del = 
                            static_cast<double>(interval_deletion_total[c][i]) /
                            interval_cov_total[c][i] / mean_deletions;
                        double interval_rel_ins = 
                            static_cast<double>(interval_insertion_total[c][i]) 
                            / interval_cov_total[c][i] / mean_insertions;
                        double interval_rel_clp = 
                            static_cast<double>(interval_clip_total[c][i]) /
                            interval_cov_total[c][i] / mean_clips;
                        double interval_rel_skp = 
                            static_cast<double>(interval_skip_total[c][i]) /
                            interval_cov_total[c][i] / mean_skips;
                        double interval_rel_start = 
                            static_cast<double>(interval_start_total[c][i]) /
                            interval_N[c][i] / mean_starts;
            
                        out << motif_number << "\t"
                            << contig_names[c] << ":" << interval_start << "-" 
                                << (interval_end + 1) << "\t"
                            << (static_cast<double>(interval_N[c][i]) /
                                genome_N * 100) << "%\t"
                            << interval_N[c][i] << "\t"
                            << 1 << "\t"
                            << interval_rel_coverage << "\t";
                        if (missing_quality)
                        {
                            out << "N/A\t";
                        }
                        else
                        {
                            out << interval_rel_quality << "\t";
                        }
                        out << interval_rel_mis << "\t"
                            << interval_rel_del << "\t"
                            << interval_rel_ins << "\t"
                            << interval_rel_clp << "\t";
                        if (missing_skips)
                        {
                            out << "N/A\t";
                        }
                        else
                        {
                            out << interval_rel_skp << "\t";
                        }
                        out << interval_rel_start << "\n";
                    }
                }
            }
            
            return;            
        }
        
        void print(ostream& out)
        {
            double mean_cov = static_cast<double>(genome_cov_total) /
                genome_N;
            double mean_qual = static_cast<double>(genome_qual_total) / 
                genome_cov_total;
            double mean_mismatch = static_cast<double>(genome_mismatch_total) / 
                genome_cov_total;
            double mean_deletions = static_cast<double>(genome_deletion_total) / 
                genome_cov_total;
            double mean_insertions = static_cast<double>(genome_insertion_total) 
                / genome_cov_total;
            double mean_clips = static_cast<double>(genome_clip_total) / 
                genome_cov_total;
            double mean_skips = static_cast<double>(genome_skip_total) /
                genome_cov_total;
            double mean_starts = static_cast<double>(genome_start_total) /
                genome_N;
   
            double motif_rel_coverage = static_cast<double>(motif_cov_total) / 
                motif_N / mean_cov;
            double motif_rel_quality = static_cast<double>(motif_qual_total) /
                motif_cov_total / mean_qual;
            double motif_rel_mis = static_cast<double>(motif_mismatch_total) /
                motif_cov_total / mean_mismatch;
            double motif_rel_del = static_cast<double>(motif_deletion_total) /
                motif_cov_total / mean_deletions;
            double motif_rel_ins = static_cast<double>(motif_insertion_total) /
                motif_cov_total / mean_insertions;
            double motif_rel_clp = static_cast<double>(motif_clip_total) /
                motif_cov_total / mean_clips;
            double motif_rel_skp = static_cast<double>(motif_skip_total) /
                motif_cov_total / mean_skips;
            double motif_rel_start = static_cast<double>(motif_start_total) /
                motif_N / mean_starts;
            
            out << std::string(80, '-');
            out << "\nMOTIF" << motif_number << " = " << motif
                << "\nMOTIF" << motif_number << " genome coverage = "
                << (static_cast<double>(motif_N) / genome_N * 100) << "%"
                << "\nMOTIF" << motif_number << " bases = " << motif_N
                << "\nMOTIF" << motif_number << " locations = "
                << total_motif_matches
                << "\nMOTIF" << motif_number << " relative coverage = "
                << motif_rel_coverage
                << "\nMOTIF" << motif_number << " relative quality = ";
            if (missing_quality)
            {
                out << "N/A";
            }
            else
            {
                out << motif_rel_quality;
            }
            out << "\nMOTIF" << motif_number << " relative mismatches = "
                << motif_rel_mis
                << "\nMOTIF" << motif_number << " relative deletions = "
                << motif_rel_del
                << "\nMOTIF" << motif_number << " relative insertions = " 
                << motif_rel_ins
                << "\nMOTIF" << motif_number << " relative clips = " 
                << motif_rel_clp
                << "\nMOTIF" << motif_number << " relative skips = ";
            if (missing_skips)
            {
                out << "N/A";
            }
            else
            {
                out << motif_rel_skp;
            }
            out << "\nMOTIF" << motif_number << " relative starts = " 
                << motif_rel_start;
            
            out << "\nMOTIF" << motif_number << " mean motif coverage = ";
            if (motif_var_size)
            {
                out << " N/A - variable size motif";
            }
            else
            {
                print_coverage_string(motif_coverage_counts,
                    mean_cov * total_motif_matches, 0, out);
            }
            out << "\n";
            
            if (mean_qual >= 0)
            {
                out << "MOTIF" << motif_number << " mean motif quality = ";
                if (motif_var_size)
                {
                    out << " N/A - variable size motif";
                }
                else
                {
                    print_quality_string(motif_quality_counts,
                        motif_coverage_counts, mean_qual, 0, out);            
                }            
                out << "\n";
            }            

            return;
        }
        
    private:
        String motif;
        size_t motif_number;
        GenomeInfoCache<bitvector>& match_cache; 
        GenomeInfoCache<MatchedLocsVec>& loc_cache;
        const vecbitvector& refamb;
        const vecbitvector& exclusions;
        const vec<String>& contig_names;
        const bool split_intervals;
        size_t current_contig;
    
        longlong motif_cov_total;
        longlong motif_qual_total;
        longlong motif_mismatch_total;
        longlong motif_deletion_total;
        longlong motif_insertion_total;
        longlong motif_clip_total;
        longlong motif_skip_total;
        longlong motif_start_total;
        longlong motif_N;
        vec<vec<longlong>> interval_cov_total;
        vec<vec<longlong>> interval_qual_total;
        vec<vec<longlong>> interval_mismatch_total;
        vec<vec<longlong>> interval_deletion_total;
        vec<vec<longlong>> interval_insertion_total;
        vec<vec<longlong>> interval_clip_total;
        vec<vec<longlong>> interval_skip_total;
        vec<vec<longlong>> interval_start_total;
        vec<vec<longlong>> interval_N;
        longlong genome_cov_total;
        longlong genome_qual_total;
        longlong genome_mismatch_total;
        longlong genome_deletion_total;
        longlong genome_insertion_total;
        longlong genome_clip_total;
        longlong genome_skip_total;
        longlong genome_start_total;
        longlong genome_N;
        longlong total_motif_matches;
        vec<longlong> motif_coverage_counts;
        vec<longlong> motif_quality_counts;
        vec<longlong> motif_start_counts;
        bool motif_var_size;
        
        bitvector* in_motif;
        MatchedLocsVec* matched_locs;
        BitVec valid_locs;
        size_t ms_loc;
        
        bool missing_quality;
        bool missing_skips;
    };
    
    /* Class for accumulating statistics for printing a coverage histogram. */
    class CoverageHistogram
    {
    public:
        CoverageHistogram(const vecbitvector& init_refamb,
            const vecbitvector& init_exclusions) :
            total_coverage(0),
            refamb(init_refamb),
            exclusions(init_exclusions)
        {
            return;
        }
        
        void process_coverage(size_t contig, size_t base,
            const CoverageAccumulator& cov)
        {
            if (!refamb[contig][base] && !exclusions[contig][base])
            {
                if (static_cast<size_t>(cov.coverage) >= coverage_count.size())
                {
                    coverage_count.resize(cov.coverage + 1, 0);
                }
                coverage_count[cov.coverage]++;
                total_coverage += cov.coverage;
                return;
            }
        }
        
        void print(ostream& out)
        {
            longlong total_locations = std::accumulate(coverage_count.begin(),
                coverage_count.end(), 0LL); // LL to invoke correct accumulate
            double mean_cov = static_cast<double>(total_coverage) / 
                total_locations;

            vec<double> bin_counts = {10, 20, 100};
            for (size_t b = 0; b < bin_counts.size(); b++)
            {
                out << std::string(80, '-') << "\n";
                out << "RELCOV\tBASES\tCUMBASES\n";
            
                longlong current_relcov_count = 0;
                longlong cum_cov_count = 0;
                for (size_t cov = 0; cov < coverage_count.size(); cov++)
                {
                    double relcov = static_cast<int>(cov / mean_cov * 
                        bin_counts[b] + 0.5) / bin_counts[b];
                
                    current_relcov_count += coverage_count[cov];
                    cum_cov_count += coverage_count[cov];
                
                    double next_relcov = static_cast<int>((cov + 1) / mean_cov
                        * bin_counts[b] + 0.5) / bin_counts[b];
                    
                    if (relcov != next_relcov ||
                        cov == coverage_count.size() - 1)
                    {
                        if (current_relcov_count > 0)
                        {
                            out << relcov << "\t"
                                << (100.0 * current_relcov_count /
                                    total_locations) << "%\t"
                                << (100.0 * cum_cov_count /
                                    total_locations) << "%\n";
                        }
                        current_relcov_count = 0;
                    }
                }
            }

            return;
        }
    
    private:
        vec<longlong> coverage_count;
        longlong total_coverage;
        const vecbitvector& refamb;
        const vecbitvector& exclusions;
    };
    
    /* Exclude any reference base that contains coverage greater
       than 99.9% of the other bases. This removes these locations from all
       bias calculations. Returns the cutoff used. */
    void exclude_coverage_outliers(vecbitvector& exclusions,
                                   BC::AbstractCoverageIterator& cov_iter,
                                   ostream& out_file)
    {
        vec<longlong> fwd_cov_hist;
        vec<longlong> rev_cov_hist;
        longlong total_fwd_counts = 0;
        longlong total_rev_counts = 0;
    
        while (cov_iter.hasNext())
        {
            BC::Coverage cov = cov_iter.next();
            if (int(fwd_cov_hist.size()) < cov.fwd + 1)
            {
                fwd_cov_hist.resize(cov.fwd + 1, 0);
            }
            if (int(rev_cov_hist.size()) < cov.rev + 1)
            {
                rev_cov_hist.resize(cov.rev + 1, 0);
            }
            fwd_cov_hist[cov.fwd]++;
            rev_cov_hist[cov.rev]++;
            total_fwd_counts++;
            total_rev_counts++;
        }
        
        longlong fwd_cum_counts = 0;
        double fwd_cum_proportion = 0;
        int fwd_cutoff = 0;
        while (fwd_cutoff < int(fwd_cov_hist.size()) &&
            fwd_cum_proportion <= 0.999)
        {
            fwd_cum_counts += fwd_cov_hist[fwd_cutoff];
            fwd_cum_proportion = double(fwd_cum_counts) / 
                double(total_fwd_counts);
            fwd_cutoff++;
        }
        fwd_cutoff += 1;
    
        longlong rev_cum_counts = 0;
        double rev_cum_proportion = 0;
        int rev_cutoff = 0;
        while (rev_cutoff < int(rev_cov_hist.size()) &&
            rev_cum_proportion <= 0.999)
        {
            rev_cum_counts += rev_cov_hist[rev_cutoff];
            rev_cum_proportion = double(rev_cum_counts) / 
                double(total_rev_counts);
            rev_cutoff++;
        }
        rev_cutoff += 1;
    
        cov_iter.rewind();
        while (cov_iter.hasNext())
        {
            BC::Coverage cov = cov_iter.next();
            if (cov.fwd > fwd_cutoff || cov.rev > rev_cutoff)
            {
                exclusions[cov.chr].Set(cov.loc, True);
            }
        }
    
        out_file << "Filtering forward coverage: max allowed was " << fwd_cutoff
                 << "\n";
        out_file << "Filtering reverse coverage: max allowed was " << rev_cutoff
                 << "\n";
    
        return;
    }
    
    /* Exclude a set of contigs. */
    bool exclude_contigs(vecbitvector& exclusions,
        const BC::GenomeReference& gen_ref, const String& ex_contigs,
        ostream& out)
    {
        vec<String> ex_contigs_vec;
        ParseStringSet(ex_contigs, ex_contigs_vec);
        
        const vec<String>& names = gen_ref.getNames();
        vec<bool> excluded_names(names.size(), false);
        size_t num_included = names.size();
        out << "excluded contigs = ";
        for (size_t e = 0; e < ex_contigs_vec.size(); e++)
        {
            bool found = false;
            DeleteLeadingWhiteSpace(ex_contigs_vec[e]);
            DeleteTrailingWhiteSpace(ex_contigs_vec[e]);
            for (size_t n = 0; !found && n < names.size(); n++)
            {
                if (names[n] == ex_contigs_vec[e])
                {
                    found = true;
                    if (!excluded_names[n])
                    {
                        exclusions[n].Set(0, exclusions[n].size(), 1);
                        if (e > 0)
                        {
                            out << ", ";
                        }
                        out << names[n];
                        excluded_names[n] = true;
                        num_included--;
                    }
                }
            }
            
            if (!found)
            {
                std::cerr << "contig " << ex_contigs_vec[e]
                    << " is not in the reference and could not be excluded\n";
                std::cerr << "available contigs:\n";
                for (size_t n = 0; n < names.size(); n++)
                {
                    std::cerr << names[n] << "\n";
                }
                return false;
            }
        }
        out << "\nincluded contigs = ";
        for (size_t n = 0; n < names.size(); n++)
        {
            if (!excluded_names[n])
            {
                out << names[n];
                num_included--;
                if (num_included > 0)
                {
                    out << ", ";
                }
            }
        }
        out << "\n";
        return true;
    }
    
    /* Exclude a set of intervals. */
    void exclude_intervals(vecbitvector& exclusions,
        const BC::GenomeReference& gen_ref, const vec<size_t>& int_chr,
        const vec<size_t>& int_start, const vec<size_t>& int_stop, ostream& out)
    {
        const vec<String>& gen_names = gen_ref.getNames();

        for (size_t i = 0; i < int_chr.size(); i++)
        {
            exclusions[int_chr[i]].Set(int_start[i], int_stop[i], 1);
        }
    
        return;
    }
    
    /* Include a set of intervals (implicitly excluding the remainder). */
    void include_intervals(vecbitvector& exclusions,
        const BC::GenomeReference& gen_ref, const vec<size_t>& int_chr,
        const vec<size_t>& int_start, const vec<size_t>& int_stop, ostream& out)
    {
        for (size_t i = 0; i < exclusions.size(); i++)
        {
            exclusions[i].Set(0, exclusions[i].size(), 1);
        }
        
        const vec<String>& gen_names = gen_ref.getNames();

        for (size_t i = 0; i < int_chr.size(); i++)
        {
            exclusions[int_chr[i]].Set(int_start[i], int_stop[i], 0);
        }
    
        return;
    }
    
    /* Include a set of contigs (implicitly excluding the remainder). */
    bool include_contigs(vecbitvector& exclusions,
        const BC::GenomeReference& gen_ref, const String& in_contigs,
        ostream& out)
    {
        vec<String> in_contigs_vec;
        if (in_contigs == "HUMAN_AUTOSOMES")
        {
            ParseStringSet("{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,"
                "20,21,22}", in_contigs_vec);
        }
        else
        {
            ParseStringSet(in_contigs, in_contigs_vec);
        }
        
        const vec<String>& names = gen_ref.getNames();
        vec<bool> names_included(names.size(), false);
        size_t count_excluded = names.size();
    
        for (size_t i = 0; i < in_contigs_vec.size(); i++)
        {
            bool found = false;
            DeleteLeadingWhiteSpace(in_contigs_vec[i]);
            DeleteTrailingWhiteSpace(in_contigs_vec[i]);
            for (size_t n = 0; !found && n < names.size(); n++)
            {
                if (names[n] == in_contigs_vec[i])
                {
                    found = true;
                    if (!names_included[n])
                    {
                        names_included[n] = true;
                        count_excluded--;
                    }
                }
            }
            
            if (!found)
            {
                std::cerr << "contig " << in_contigs_vec[i]
                    << " is not in the reference and could not be included\n";
                std::cerr << "available contigs:\n";
                for (size_t n = 0; n < names.size(); n++)
                {
                    std::cerr << names[n] << "\n";
                }
                return false;
            }
        }
        
        if (count_excluded > 0)
        {
            String exclusion_set = "{";
            for (size_t n = 0; n < names.size(); n++)
            {
                if (!names_included[n])
                {
                    exclusion_set += names[n] + ",";
                }
            }
            exclusion_set[exclusion_set.size() - 1] = '}';
            exclude_contigs(exclusions, gen_ref, exclusion_set, out);
        }
    
        return true;
    }
    
    /* If the user has specified a predefined motif set, return that set of
       motifs. */
    String fetch_motif_sets(const String& set_name)
    {
        if (set_name == "STANDARD")
        {
            return "{N^50 (G+C >= 0.85)^100 N^50,"
                "N^50 (G+C >= 0.75)^100 N^50,"
                "N^50 (G >= 0.8)^30 N^50,"
                "N^50 (AT)^15 N^50,"
                "N^50 (G+C <= 0.1)^100 N^50}";
        }
        else if (set_name == "WEAKSTANDARD")
        {
            return "{N^50 (G+C >= 0.70)^100 N^50,"
                "N^50 (G >= 0.75)^30 N^50,"
                "N^55 (AT)^10 N^55,"
                "N^50 (G+C <= 0.13)^100 N^50}";
        }
        else if (set_name == "HOMOPOLYMERS")
        {
            String homopolymer_motif_set = "{";
            for (int i = 1; i <= 20; i++)
            {
                homopolymer_motif_set += "N (HOMOPOLYMER)^" + ToString(i) + 
                    " N";
                
                if (i != 20)
                {
                    homopolymer_motif_set += ",";
                }
            }
            homopolymer_motif_set += "}";
            return homopolymer_motif_set;
        }
        else if (set_name == "NULL")
        {
            return "";
        }
        else
        {
            return set_name;
        }
    }
} // end nameless namespace

/* Main function of BadCoverage, with many command-line options. */
int main(int argc, char *argv[])
{
    RunTime();

    BeginCommandArgumentsNoHeader;
    CommandArgument_String_Doc(COV, "a (possibly mixed) "
    	"space-delimited list of BAM/SCI filenames or Picard lane specifiers "
    	"(e.g. 61CCEAAXX.2 or, more specifically, to select a particular "
    	"library in an indexed lane, append that, e.g. "
    	"61CCEAAXX.2.Solexa-38535; to select a particular "
    	"run of a particular lane, append that, e.g. "
    	"61CEEAAXX.2.Solexa-38535.C1-202_2010-09-15_2010-10-18 or "
    	"61CEEAAXX.2..C1-202_2010-09-15_2010-10-18 if you don't wish to "
    	"specify the library manually). Specify an empty string to just "
    	"populate the cache.");
    CommandArgument_String_OrDefault_Doc(CACHE_DIR,
        getenv("HOME") + String("/.badcoverage/cache"),
        "to cache vecbitvector of bases marked by motif on genome");
    CommandArgument_String_OrDefault_Doc(REF_FASTA, "",
        "fasta file for reference, if absent, will try to infer from BAM "
        "header");
    CommandArgument_DoubleSet_OrDefault_Doc(FRAC_TO_USE, "1",
        "randomly discard aligned reads to get to this fraction of total");
    CommandArgument_Bool_OrDefault_Doc(RAND_CLOCK_SEED, True,
        "seed the random number generator from the clock so runs with the "
        "same FRAC_TO_USE percentage aren't identical");
    CommandArgument_String_OrDefault_Doc(TEST_INTERVALS, "",
        "file of intervals, one line per interval, form is chr*:*-*; "
        "compare coverage on these intervals to overall coverage; "
        "multiple files can be specified with set notation");
    CommandArgument_String_OrDefault_Doc(TEST_MOTIFS, "STANDARD",
        "list of motifs - compare coverage on these intervals to overall "
        "coverage");
    CommandArgument_String_OrDefault_Doc(COVERAGE_MOTIFS_OUT, "",
        "prefix of filename to write coverage motif intervals");
    CommandArgument_String_OrDefault_Doc(MOTIF_OUT, "", "prefix of files "
        "that will contain locations");
    CommandArgument_Double_OrDefault_Doc(BAD_COVERAGE_THRESH, 0.1,
        "relative coverage less than this percentage indicates a bad "
        "location");
    CommandArgument_Int_OrDefault_Doc(EXCLUDE_MAPQ_LTE, -1,
        "exclude reads <= this mapping quality threshold");
    CommandArgument_String_OrDefault_Doc(OUTPUT_FILE, "/dev/stdout",
        "file to direct non-error output to, defaults to standard output");
    CommandArgument_Bool_OrDefault_Doc(IGNORE_DUPLICATES, True,
        "ignore reads marked as duplicates");
    CommandArgument_Bool_OrDefault_Doc(IGNORE_NONPF, True,
        "ignore reads that fail purity filtering");
    CommandArgument_Int_OrDefault_Doc(COVERAGE_CAP, 0,
        "if non-zero, the maximum amount of per-base coverage counted");
    CommandArgument_Bool_OrDefault_Doc(ONLY_READ1, False,
        "Only count read 1 in pairs");
    CommandArgument_Bool_OrDefault_Doc(ONLY_READ2, False,
        "Only count read 2 in pairs");
    CommandArgument_Bool_OrDefault_Doc(ROBUST_COVERAGE, False,
        "Exclude regions with coverage greater than 99.9% of genome "
        "locations from the calculation");
    CommandArgument_Int_OrDefault_Doc(PAIR_FILTERING, 0,
        "if non-zero, any pairs not aligned within this distance of one "
        "another will be ignored");
    CommandArgument_Bool_OrDefault_Doc(LATENT_VAR_OUTPUT, False,
        "Output an estimate of the latent variance of read-start probabilities "
        "across the reference");
    CommandArgument_Bool_OrDefault_Doc(LATENT_VAR_CI_OUTPUT, False,
        "Output an (currently too large) 95% confidence interval for "
        "the latent variable statistic");
    CommandArgument_String_OrDefault_Doc(EXCLUDE_CONTIGS, "", "Exclude this "
        "set of contigs from the statistics, even if they include intervals "
        "specified via the INCLUDE_INTERVALS argument.");
    CommandArgument_String_OrDefault_Doc(INCLUDE_CONTIGS, "", "Exclude the " 
        "contigs that do *not* appear in this set.");
    CommandArgument_String_OrDefault_Doc(EXCLUDE_INTERVALS, "", "Exclude the "
        "intervals listed in this file from the statistics.");
    CommandArgument_String_OrDefault_Doc(INCLUDE_INTERVALS, "", "Include the "
        "intervals listed in this file from the statistics, unless they "
        "conflict with EXCLUDE_CONTIGS or INCLUDE_CONTIGS.");
    CommandArgument_Bool_OrDefault_Doc(ENABLE_COVERAGE_MOTIFS, False,
        "Compute bad fraction and other coverage motifs, only suitable for "
        "deep (~100x) coverage data.")
    CommandArgument_UnsignedInt_OrDefault_Doc(GC_BIAS_TABLE_WINSIZE, 100,
        "Print a tab-delimited table of relative coverage, quality, and error "
        "values for non-overlapping GC_BIAS_TABLE_WINSIZE-base GC windows. "
        "The table can be disabled by setting GC_BIAS_TABLE_WINSIZE=0.");
    CommandArgument_Bool_OrDefault_Doc(PRINT_COVERAGE_HISTOGRAM, False,
        "Print a histogram indicating what percent of bases are at each "
        "relative coverage level.");
    CommandArgument_String_OrDefault_Doc(SERIALIZE_COVERAGE, "", "Save the "
        "coverage data to a binary file to enable fast (high-memory) "
        "reprocessing.");
    CommandArgument_Int_OrDefault_Doc(INSERT_TEST_INTERVALS_INDEX, -1,
        "MOTIF number of first test interval (with the remainder to be "
        "inserted in sequence afterwards) - the default is for test "
        "intervals to be added after the final TEST_MOTIF.");
    CommandArgument_Bool_OrDefault_Doc(USE_ORIGINAL_QSCORES, True,
        "Use the original quality scores if present in a SAM/BAM input.");
    CommandArgument_Bool_OrDefault_Doc(USE_TABLE_FORMAT, False,
        "Format the motif reports as a spreadsheet, tab-separated table.");
    CommandArgument_Bool_OrDefault_Doc(SPLIT_MOTIFS, False,
        "Produce sub-reports on each matching interval of each motif, probably "
        "best combined with USE_TABLE_FORMAT.");
    EndCommandArguments;

    std::ofstream nonstd_ofstream;
    
    if (OUTPUT_FILE != "/dev/stdout")
    {
        nonstd_ofstream.open(OUTPUT_FILE.c_str());
    }
    
    std::ostream& out_file = (OUTPUT_FILE == "/dev/stdout") ? std::cout : 
        nonstd_ofstream;

    command.PrintTheCommandPretty(out_file);

    // Check arguments.
    for (size_t f = 0; f < FRAC_TO_USE.size(); f++)
    {
        ForceAssertGe(FRAC_TO_USE[f], 0);
        ForceAssertLe(FRAC_TO_USE[f], 1);
    }
    
    if (SPLIT_MOTIFS && !USE_TABLE_FORMAT)
    {
        FatalErr("SPLIT_MOTIFS=True must be combined with "
            " USE_TABLE_FORMAT=True");
    }
    
    // discover the BAM/SCI files to process
    vec<String> cov_files;
    if (!COV.empty())
    {
        cov_files = BC::locate_cov_files(COV);
        out_file << "coverage files surveyed = ";
        CompactPrint(out_file, cov_files, ", ");
        out_file << "\n";
    }
    
    if (FRAC_TO_USE.size() != 1 && FRAC_TO_USE.size() != cov_files.size())
    {
        FatalErr("FRAC_TO_USE must contain one value or a list equal to the "
            "number of coverage files");
    }
    
    // find reference and load the relevant info
    if (REF_FASTA.empty())
    {
        if (cov_files.size() > 0)
        {
            REF_FASTA = BC::locate_reference(cov_files);

            if (REF_FASTA.empty())
            {
                std::cerr << "Failed to locate reference file.\n";
                return EXIT_FAILURE;
            }
        }
        else
        {
            std::cerr << "You need to specify a REF_FASTA.\n";
            return EXIT_FAILURE;
        }
    }
    
    if (!IsRegularFile(REF_FASTA))
    {
        std::cerr << REF_FASTA << " is not a regular file.\n";
        return EXIT_FAILURE;
    }
    else
    {
        // canonicalize the reference FASTA's path, to help us find the right
        // cache directory
        REF_FASTA = RealPath(REF_FASTA);
        out_file << "reference = " << REF_FASTA << "\n";
    }
    
    CACHE_DIR = make_cache_dir(CACHE_DIR, REF_FASTA, out_file);

    BC::GenomeReference gen_ref(REF_FASTA, CACHE_DIR);
    const vecbitvector& refamb = gen_ref.getAmbiguities();
    const vecbasevector& genome = gen_ref.getBases();
    const vec<String>& names = gen_ref.getNames();
    
    // Create motifs, if necessary
    String motifmatch_cache_dir = CACHE_DIR + "/MOTIFMATCH";
    String motifloc_cache_dir = CACHE_DIR + "/MOTIFLOC";

    TEST_MOTIFS = fetch_motif_sets(TEST_MOTIFS);
    
    vec<String> motifs = parse_motifs(TEST_MOTIFS);
    vec<String> motif_dir_names = initialize_motif_caches
        (motifmatch_cache_dir, motifloc_cache_dir, motifs, gen_ref, out_file);

    vec<GenomeInfoCache<bitvector>* > match_caches(motifs.size());
    vec<GenomeInfoCache<MatchedLocsVec>* > loc_caches(motifs.size());
    
    for (size_t mi = 0; mi < match_caches.size(); mi++)
    {
        match_caches[mi] = new GenomeInfoCache<bitvector>(genome,
            motifmatch_cache_dir + "/" + motif_dir_names[mi], false, true);
        loc_caches[mi] = new GenomeInfoCache<MatchedLocsVec>(genome,
            motifloc_cache_dir + "/" + motif_dir_names[mi], false, true);
    }

    // Load test intervals.
    if (TEST_INTERVALS.empty() && !CACHE_DIR.empty())
    {
        String cached_ti_file = CACHE_DIR + "/BadCoverage_test_intervals";
        if (IsRegularFile(cached_ti_file))
        {
            TEST_INTERVALS = cached_ti_file;
        }
    }
    
    if (!TEST_INTERVALS.empty() && TEST_INTERVALS != "NULL")
    {
        vec<String> ti_files;
        if (TEST_INTERVALS.StartsWith("{") && TEST_INTERVALS.EndsWith("}"))
        {
            ParseStringSet(TEST_INTERVALS, ti_files);
        }
        else
        {
            ti_files.push_back(TEST_INTERVALS);
        }

        size_t insert_intervals_base = motifs.size();
        if (INSERT_TEST_INTERVALS_INDEX >= 0)
        {
            insert_intervals_base = INSERT_TEST_INTERVALS_INDEX;
        }
        
        if (insert_intervals_base > motifs.size())
        {
            std::cerr << "INSERT_TEST_INTERVALS_INDEX cannot exceed the "
                "number of motifs\n";
            return EXIT_FAILURE;
        }
        
        for (size_t t = 0; t < ti_files.size(); t++)
        {
            DeleteLeadingWhiteSpace(ti_files[t]);
            DeleteTrailingWhiteSpace(ti_files[t]);
            
            bool inverted = false;
            if (ti_files[t][0] == '-')
            {
                inverted = true;
                ti_files[t] = ti_files[t].substr(1, ti_files[t].size());
            }
            
            out_file << "Test Intervals [" << t << "] file = " << ti_files[t]
                << (inverted ? " (Inverted)\n" : "\n");
            
            vec<size_t> test_chr(0);
            vec<size_t> test_start(0);
            vec<size_t> test_stop(0);
            load_intervals(test_chr, test_start, test_stop, gen_ref,
                ti_files[t]);
            
            match_caches.insert(match_caches.begin() + insert_intervals_base + 
                t, 0);
            loc_caches.insert(loc_caches.begin() + insert_intervals_base + t, 
                0);
            motifs.insert(motifs.begin() + insert_intervals_base + t, 
                String("Test Intervals [" + ToString(t) + "]"
                + (inverted ? " (Inverted)" : "")));
            
            if (inverted)
            {
                testintervals2inversemotif(test_chr, test_start, test_stop, 
                    gen_ref, match_caches[insert_intervals_base + t],    
                    loc_caches[insert_intervals_base + t]);            
            }
            else
            {
                testintervals2motif(test_chr, test_start, test_stop, gen_ref,
                    match_caches[insert_intervals_base + t],    
                    loc_caches[insert_intervals_base + t]);
            }
        }
    }
    
    if (MOTIF_OUT != "")
    {
        for (size_t mi = 0; mi < motifs.size(); mi++)
        {
            Ofstream(mout, MOTIF_OUT + ToString(mi) + ".intervals");
            for (size_t t = 0; t < genome.size(); t++)
            {
                bitvector* matchvec = match_caches[mi]->getInfo(t);
                for (size_t j = 0; j < matchvec->size(); j++)
                {
                    if ((*matchvec)[j])
                    {
                        size_t k = matchvec->NextDiff(j);
                        mout << names[t] << ":" << j << "-" << k << "\n";
                        j = k - 1;
                    }
                }
            }
            mout.close();
        }
    }
    
    // quit if there are no COVs - the user just wanted to create (or add
    // motifs to) the persistent motif cache for the reference
    if (cov_files.size() == 0)
    {
        out_file << "Motif caches created. No COV argument was specified, "
                 << "so BadCoverage has completed its assignment.\n";
        return EXIT_SUCCESS;
    }
    
    BC::AbstractCoverageIterator* aciter =
        new BC::MultiCoverageIterator(cov_files, gen_ref,
        IGNORE_DUPLICATES, IGNORE_NONPF, EXCLUDE_MAPQ_LTE, FRAC_TO_USE,
        RAND_CLOCK_SEED, PAIR_FILTERING, ONLY_READ1, ONLY_READ2, COVERAGE_CAP,
        USE_ORIGINAL_QSCORES);

    vecbitvector exclusions;
    Mimic(genome, exclusions);
    // first mark up specific intervals    
    if (!INCLUDE_INTERVALS.empty() && !EXCLUDE_INTERVALS.empty())
    {
        std::cerr << "Only at most one of INCLUDE_INTERVALS and "
            "EXCLUDE_INTERVALS can be specified.\n";
        return EXIT_FAILURE;
    }
    else if (!INCLUDE_INTERVALS.empty())
    {
        vec<size_t> in_int_chr(0);
        vec<size_t> in_int_start(0);
        vec<size_t> in_int_stop(0);
        load_intervals(in_int_chr, in_int_start, in_int_stop, gen_ref,
            INCLUDE_INTERVALS);
        include_intervals(exclusions, gen_ref, in_int_chr, in_int_start,
            in_int_stop, out_file);
        out_file << "included intervals = " << INCLUDE_INTERVALS << "\n";
    }
    else if (!EXCLUDE_INTERVALS.empty())
    {
        vec<size_t> ex_int_chr(0);
        vec<size_t> ex_int_start(0);
        vec<size_t> ex_int_stop(0);
        load_intervals(ex_int_chr, ex_int_start, ex_int_stop, gen_ref,
            EXCLUDE_INTERVALS);
        exclude_intervals(exclusions, gen_ref, ex_int_chr, ex_int_start,
            ex_int_stop, out_file);
        out_file << "excluded intervals = " << EXCLUDE_INTERVALS << "\n";
    }
    // now exclude excluded (or non-included) contigs
    if (!INCLUDE_CONTIGS.empty() && !EXCLUDE_CONTIGS.empty())
    {
        std::cerr << "Only at most one of INCLUDE_CONTIGS and EXCLUDE_CONTIGS "
            "can be specified.\n";
        return EXIT_FAILURE;
    }    
    else if (!EXCLUDE_CONTIGS.empty())
    {
        if (!exclude_contigs(exclusions, gen_ref, EXCLUDE_CONTIGS, out_file))
        {
            std::cerr << "Excluding contigs failed.\n";
            return EXIT_FAILURE;
        }
    }
    else if (!INCLUDE_CONTIGS.empty())
    {
        if (!include_contigs(exclusions, gen_ref, INCLUDE_CONTIGS, out_file))
        {
            std::cerr << "Including contigs failed.\n";
            return EXIT_FAILURE;
        }
    }
    
    // information about the libraries, production units (e.g. flowcells)
    // and samples involved    
    out_file << "libraries = ";
    CompactPrint(out_file, aciter->getLibraries(), ", ");
    out_file << "\n";
    
    out_file << "production units = ";
    CompactPrint(out_file, aciter->getProductionUnits(), ", ");
    out_file << "\n";

    out_file << "samples = ";
    CompactPrint(out_file, aciter->getSamples(), ", ");
    out_file << "\n";
        
    if (ROBUST_COVERAGE)
    {
        // run through the BAM files and find the unweighted coverage
        // distribution - this is then used to exclude genome
        // locations with abnormally high coverage (they are marked
        // as "ambiguous")
        exclude_coverage_outliers(exclusions, *aciter, out_file);
        aciter->rewind();
    }

    BC::SavedCoverageIterator* sciter = NULL;
    
    if (!SERIALIZE_COVERAGE.empty())
    {
        sciter = new BC::SavedCoverageIterator(SERIALIZE_COVERAGE,
            REF_FASTA, aciter->getLibraries(), aciter->getProductionUnits(),
            aciter->getSamples());
    }

    // analyze data
    CoverageAccumulator covacc;
    BulkStatistics bulkstats(refamb, exclusions);
    LatentVarianceAnalysis lvanalysis(LATENT_VAR_OUTPUT, 
        LATENT_VAR_CI_OUTPUT, refamb, exclusions);
    vec<ContentMotifStatistics*> motifstats(motifs.size());
    for (size_t mi = 0; mi < motifstats.size(); mi++)
    {
        motifstats[mi] = new ContentMotifStatistics(motifs[mi], mi,
            *match_caches[mi], *loc_caches[mi], refamb, exclusions, names,
            SPLIT_MOTIFS);
    }
    GCBiasTable gcbtable(GC_BIAS_TABLE_WINSIZE, genome, refamb, exclusions);
    CoverageHistogram covhist(refamb, exclusions);
            
    for (size_t t = 0; t < refamb.size(); t++)
    {
        for (size_t i = 0; i < refamb[t].size(); i++)
        {
            if (!aciter->hasNext())
            {
                std::cerr << "ERROR: CoverageIterator ended before "
                    << "reference at " << names[t] << ":" << i
                    << ".\n";
                CRD::exit(EXIT_FAILURE);
            }

            BC::Coverage cov = aciter->next();

            if (!SERIALIZE_COVERAGE.empty())
            {
                sciter->addCoverage(cov, aciter->numReads(),
                    aciter->numReadsUsed());
            }
            
            if (cov.chr != t || cov.loc != i)
            {
                std::cerr << "cov.chr = " << cov.chr << " t = " << t
                    << "\n";
                std::cerr << "cov.loc = " << cov.loc << " i = " << i
                    << "\n";
                std::cerr << "ERROR: CoverageIterator not synchronized "
                    << "with reference.\n";
                CRD::exit(EXIT_FAILURE);
            }
            
            covacc += cov;

            if (!ENABLE_COVERAGE_MOTIFS) // otherwise, compute on second pass
            {
                for (size_t mi = 0; mi < motifstats.size(); mi++)
                {
                    motifstats[mi]->process_coverage(t, i, covacc);
                }
            }

            bulkstats.process_coverage(t, i, covacc);
            gcbtable.process_coverage(t, i, covacc);
            lvanalysis.process_coverage(t, i, covacc);
            covhist.process_coverage(t, i, covacc);
        }
    }

    if (aciter->hasNext())
    {
        std::cerr << "ERROR: Reference ended before CoverageIterator\n";
        CRD::exit(EXIT_FAILURE);
    }

    if (!SERIALIZE_COVERAGE.empty())
    {
        delete sciter;
    }

    vec<CoverageMotif*> covmotifs;
    if (ENABLE_COVERAGE_MOTIFS)
    {
        if (!SERIALIZE_COVERAGE.empty())
        {
            // cleverly we will switch to using our recently serialized
            // coverage for the next pass, probably saving a lot of time
            delete aciter;
            aciter = new BC::SavedCoverageIterator(SERIALIZE_COVERAGE);
        }
        else
        {
            aciter->rewind();
        }
        
        covmotifs.push_back(new BadCoverageMotif(BAD_COVERAGE_THRESH * 
            bulkstats.mean_coverage(), motifs, refamb, exclusions));
        covmotifs.push_back(new InverseMotif(*covmotifs[0],
            "good coverage", motifs, refamb, exclusions));
        if (!bulkstats.missing_quality())
        {
            covmotifs.push_back(new BadQualityMotif(bulkstats.mean_quality()
                + 10 * log10(BAD_COVERAGE_THRESH), motifs, refamb, exclusions));
        }
        covmotifs.push_back(new IntersectionMotif(*covmotifs.back(), 
            *covmotifs[1], motifs, refamb, exclusions));
        covmotifs.push_back(new BadMismatchesMotif(bulkstats.mean_mismatches()
            / BAD_COVERAGE_THRESH, motifs, refamb, exclusions));
        covmotifs.push_back(new IntersectionMotif(*covmotifs.back(), 
            *covmotifs[1], motifs, refamb, exclusions));
        covmotifs.push_back(new BadDeletionsMotif(bulkstats.mean_deletions()
            / BAD_COVERAGE_THRESH, motifs, refamb, exclusions));
        covmotifs.push_back(new IntersectionMotif(*covmotifs.back(), 
            *covmotifs[1], motifs, refamb, exclusions));
        covmotifs.push_back(new BadInsertionsMotif(bulkstats.mean_insertions()
            / BAD_COVERAGE_THRESH, motifs, refamb, exclusions));
        covmotifs.push_back(new IntersectionMotif(*covmotifs.back(), 
            *covmotifs[1], motifs, refamb, exclusions));
        covmotifs.push_back(new BadClipsMotif(bulkstats.mean_clips() / 
            BAD_COVERAGE_THRESH, motifs, refamb, exclusions));
        covmotifs.push_back(new IntersectionMotif(*covmotifs.back(), 
            *covmotifs[1], motifs, refamb, exclusions));
        covmotifs.push_back(new UncoveredMotif(motifs, refamb, exclusions));
        
        vec<bool> genome_motif_matched(motifs.size(), false);
        for (size_t t = 0; t < refamb.size(); t++)
        {
            for (size_t i = 0; i < refamb[t].size(); i++)
            {
                if (!aciter->hasNext())
                {
                    std::cerr << "ERROR: CoverageIterator ended before "
                        << "reference at " << names[t] << ":" << i
                        << ".\n";
                    CRD::exit(EXIT_FAILURE);
                }
    
                BC::Coverage cov = aciter->next();
                
                if (cov.chr != t || cov.loc != i)
                {
                    std::cerr << "cov.chr = " << cov.chr << " t = " << t
                        << "\n";
                    std::cerr << "cov.loc = " << cov.loc << " i = " << i
                        << "\n";
                    std::cerr << "ERROR: CoverageIterator not synchronized "
                        << "with reference.\n";
                    CRD::exit(EXIT_FAILURE);
                }
                
                covacc += cov;
                    
                for (size_t mi = 0; mi < motifstats.size(); mi++)
                {
                    genome_motif_matched[mi] =
                        motifstats[mi]->process_coverage(t, i, covacc);
                }
                
                for (size_t cm = 0; cm < covmotifs.size(); cm++)
                {
                    covmotifs[cm]->process_coverage(t, i, covacc,
                        genome_motif_matched);                            
                }
            }
        }
    }

    longlong num_ref_bases = 0;
    for (size_t g = 0; g < genome.size(); g++)
    {
        num_ref_bases += genome[g].size();
    }
        
    // print report
    out_file << "total number of reads/pairs = " << aciter->numReads() << "\n";
    out_file << "number of mapped reads = " << aciter->numReadsUsed() << "\n";
    out_file << "reads per base = " << aciter->numReadsUsed() /
        static_cast<double>(num_ref_bases) << "\n";    
    
    bulkstats.print(out_file);
    lvanalysis.print(out_file);
        
    for (size_t cm = 0; cm < covmotifs.size(); cm++)
    {
        covmotifs[cm]->print(out_file);
        if (COVERAGE_MOTIFS_OUT != "")
        {
            String coverage_motifs_out_fname = FilenameSafeString
                (COVERAGE_MOTIFS_OUT + "_"
                + covmotifs[cm]->name() + ".intervals");

            Ofstream(cmout, coverage_motifs_out_fname);
            covmotifs[cm]->print_intervals(cmout, names);
            cmout.close();
            
            String coverage_motifs_nm_out_fname = FilenameSafeString
                (COVERAGE_MOTIFS_OUT + "_"
                + covmotifs[cm]->name() + "_nonmotif.intervals");
            Ofstream(cmnonmotifout, coverage_motifs_nm_out_fname);
            covmotifs[cm]->print_nongenomemotif_intervals
                (cmnonmotifout, names);
            cmnonmotifout.close();
        }
        
        delete covmotifs[cm];
    }    
    
    if (USE_TABLE_FORMAT)
    {
        ContentMotifStatistics::print_table_header(out_file);
    }
    
    for (size_t mi = 0; mi < motifs.size(); mi++)
    {
        if (USE_TABLE_FORMAT)
        {
            motifstats[mi]->print_table_data(out_file);
        }
        else
        {
            motifstats[mi]->print(out_file);
        }
        delete motifstats[mi];
    }

    if (GC_BIAS_TABLE_WINSIZE)
    {
        gcbtable.print(out_file);
    }
    
    if (PRINT_COVERAGE_HISTOGRAM)
    {
        covhist.print(out_file);
    }

    delete aciter;

    for (size_t mi = 0; mi < match_caches.size(); mi++)
    {
        delete match_caches[mi];
        delete loc_caches[mi];
    }

    out_file.flush();

    if (nonstd_ofstream.is_open())
    {
        nonstd_ofstream.close();
    }
    
    return EXIT_SUCCESS;
}
