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
#include "Basevector.h"
#include "Bitvector.h"
#include "FetchReads.h"
#include "FetchReadsAmb.h"
#include "bias/CoverageIterator.h"
#include "bias/GenomeReference.h"
#include "system/AllOfOutput.h"

namespace BC
{
    /* Load genome information from a reference FASTA. */
    GenomeReference::GenomeReference(const String& init_ref_fasta,
        const String& cache_dir) : ref_fasta(init_ref_fasta)
    {
        String names_cache = cache_dir + "/BadCoverage_names";
        String refamb_cache = cache_dir + "/BadCoverage_refamb";
        String genome_cache = cache_dir + "/BadCoverage_genome";
        
        if (cache_dir != "" && IsRegularFile(names_cache)
            && IsRegularFile(refamb_cache) && IsRegularFile(genome_cache))
        {
            if (!IsRegularFile(cache_dir + "/finished.txt"))
            {
                std::cerr << "Reference cache " << cache_dir << 
                    " not completely initialized - either some "
                    "other process is building it or it is "
                    "corrupt and should be deleted.\n";
                CRD::exit(EXIT_FAILURE);
            }
            BinaryReader::readFile(names_cache.c_str(),&names);
            ambiguities.ReadAll(refamb_cache);
            bases.ReadAll(genome_cache);
        }
        else
        {
            vecString vs_names;
            FetchReads(bases, vs_names, ref_fasta);
            names.insert(names.begin(), vs_names.begin(), vs_names.end());
            for (size_t n = 0; n < names.size(); n++)
            {
                if (names[n].Contains(" "))
                {
                    names[n] = names[n].Before(" ");
                }
            }
            FetchReadsAmb(ambiguities, ref_fasta);
    
            if (cache_dir != "")
            {
                if (!IsRegularFile(names_cache)
                    && !IsRegularFile(refamb_cache)
                    && !IsRegularFile(genome_cache))
                {
                    if (!IsRegularFile(cache_dir + "/working.txt"))
                    {
                        std::ofstream working((cache_dir +
                            "/working.txt").c_str());
                        working << "caching " << ref_fasta << "\n";
                        working.close();
                    }
                    else
                    {
                        std::cerr << "Another process is building a "
                            "reference cache in " << cache_dir << "\n";
                        CRD::exit(1);
                    }
                    
                    BinaryWriter::writeFile(names_cache.c_str(),names);
                    ambiguities.WriteAll(refamb_cache);
                    bases.WriteAll(genome_cache);
                    
                    Remove(cache_dir + "/working.txt");
                    std::ofstream finished((cache_dir +
                        "/finished.txt").c_str());
                    finished << ref_fasta << " cached\n";
                    finished.close();
                }
                else
                {
                    std::cerr << "CACHE_DIR " << cache_dir
                        << " is inconsistent - there should be "
                        << names_cache << ", " << refamb_cache << ", and "
                        << genome_cache << " files, or none of them."
                        << std::endl;
                    CRD::exit(EXIT_FAILURE);
                }
            }
        }
    
        return;    
    }
    
    /* Access the bases (N's are represented by a random base). */
    const vecbasevector& GenomeReference::getBases() const
    {
        return bases;
    }
    
    /* Access bitvectors indicating ambiguous bases. */
    const vecbitvector& GenomeReference::getAmbiguities() const
    {
        return ambiguities;
    }
    
    /* Access contig names. */
    const vec<String>& GenomeReference::getNames() const
    {
        return names;
    }
    
    /* Look up a reference file using a Picard directory's params.txt files
       The argument is a vector of Picard-located BAMs and the function returns
       an empty string on failure. */
    String locate_reference_from_lane(const String& picard_bam)
    {
        String picard_dir = picard_bam.RevBefore("/");
        vec<String> fastas = AllOfOutput("grep '^REFERENCE_SEQUENCE=' "
            + picard_dir + "/params.txt");
        if (fastas.size() != 1)
        {
            std::cerr << picard_dir << " has a misformatted params.txt file"
                << std::endl;
            CRD::exit(EXIT_FAILURE);
        }
        else
        {
            return fastas[0].After("REFERENCE_SEQUENCE=");
        }
    }
    
    /* Look up a reference file using a SCI file header. */
    String locate_reference_from_sci(const String& sci_filename)
    {
        BC::SavedCoverageIterator sci(sci_filename);
        return sci.getReference();
    }
    
    /* Look up a reference file using a BAM file header. */
    String locate_reference_from_bam(const String& bam_filename)
    {
        SAM::BAMFile sf(bam_filename,"",true);
        String refloc = sf.getUniqueRefURI();
    
        if (refloc.Contains("file:", 0))
        {
            refloc = refloc.After("file:");
        }
        
        return refloc;
    }
    
    /* Look up a reference file using BAM/SCI file headers - check that all 
       refer to the same FASTA. Return empty string on failure. */
    String locate_reference(const vec<String>& cov_files)
    {
        vec<String> ref_fastas;
        for (size_t j = 0; j < cov_files.size(); j++)
        {
            if (cov_files[j].StartsWith("/seq/") && 
                cov_files[j].Contains("/picard/"))
            {
                ref_fastas.push_back(locate_reference_from_lane(cov_files[j]));
            }
            else if (cov_files[j].EndsWith(".bam") || 
                cov_files[j].EndsWith(".sam"))
            {
                ref_fastas.push_back(locate_reference_from_bam(cov_files[j]));
            }
            else if (cov_files[j].EndsWith(".sci"))
            {
                ref_fastas.push_back(locate_reference_from_sci(cov_files[j]));
            }
            else
            {
                std::cerr << "Unknown coverage file format: " + cov_files[j]
                    << "\n";
                return "";
            }
        }
        UniqueSort(ref_fastas);
        if (!ref_fastas.solo())
        {
            std::cerr << "Not sure which reference file to use.\n";
            std::cerr << "These reference files were found:\n";
            for (size_t f = 0; f < ref_fastas.size(); f++)
            {
                std::cerr << ref_fastas[f] << std::endl;
            }
            return "";
        }
    
        return ref_fastas[0];
    }
}
