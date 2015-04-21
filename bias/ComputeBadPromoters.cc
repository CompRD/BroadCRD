///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// This program computes the "bad promoter" list given coverage data and the
// RefSeq database. It is designed to replicate the method described by Gad
// Getz and Michael Lawrence and (hopefully) published in the methods section
// of "Characterizing and measuring bias in sequence data" by Ross et al.

#include "FeudalMimic.h"
#include "Intvector.h"
#include "MainTools.h"
#include "VecUtilities.h"
#include "bias/CoverageIterator.h"
#include "bias/GenomeReference.h"

int main(int argc, char** argv)
{
    RunTime();
    
    BeginCommandArguments;
    CommandArgument_String_Doc(COV, "coverage BAM/SCI file");
    CommandArgument_String_Doc(REF, "genome reference file");
    CommandArgument_String_Doc(REFSEQDB, "REFSEQ database for transcription "
        "start sites");
    CommandArgument_String_OrDefault_Doc(REF_CACHE_DIR,
        getenv("HOME") + String("/.badcoverage/cache"),
        "directory of cached genome reference information");
    EndCommandArguments;

    REF = RealPath(REF);
    if (!REF_CACHE_DIR.empty())
    {
        REF_CACHE_DIR = REF_CACHE_DIR + "/" + REF;
    }
    
    BC::GenomeReference gref(REF, REF_CACHE_DIR);
    std::cout << "reference loaded" << std::endl;
    BC::AbstractCoverageIterator* coviter = new 
        BC::MultiCoverageIterator(BC::locate_cov_files(COV), gref);
    std::cout << "coverage iterator initialized" << std::endl;

    vec<String> names = gref.getNames();
    const vecbitvector& refamb = gref.getAmbiguities();
    const vecbasevector& genome = gref.getBases();
    
    vec<size_t> name_ids(names.size());
    for (size_t n = 0; n < name_ids.size(); n++)
    {
        name_ids[n] = n;
    }
    SortSync(names, name_ids);
    
    VecIntVec coverage;
    Mimic(genome, coverage);

    std::cout << "building coverage map" << std::endl;
    while (coviter->hasNext())
    {
        BC::Coverage c = coviter->next();
        coverage[c.chr][c.loc] = c.fwd + c.rev;
    }
    
    std::ifstream refseq_stream(REFSEQDB);
    String refseq_header;
    getline(refseq_stream, refseq_header);

    // parse RefSeq header for necessary field locations
    vec<String> headers;
    Tokenize(refseq_header, '\t', headers);
    
    int chrom = -1;
    int txStart = -1;
    int txEnd = -1;
    int strand = -1;
    int name2 = -1;
    
    for (size_t h = 0; h < headers.size(); h++)
    {
        if (headers[h] == "chrom")
        {
            chrom = h;
        }
        else if (headers[h] == "txStart")
        {
            txStart = h;
        }
        else if (headers[h] == "txEnd")
        {
            txEnd = h;
        }
        else if (headers[h] == "strand")
        {
            strand = h;
        }
        else if (headers[h] == "name2")
        {
            name2 = h;
        }
    }
    
    if (chrom == -1 || txStart == -1 || txEnd == -1 || strand == -1 ||
        name2 == -1)
    {
        FatalErr("Unable to parse RefSeq header");
    }
    vec<String> bad_promoter_intervals;
    vec<double> bad_promoter_ratios;
    vec<String> bad_promoter_genes;
    while (refseq_stream)
    {
        String refseq_record;
        getline(refseq_stream, refseq_record);
        if (refseq_record.size() > 0)
        {
            vec<String> refseq_fields;
            Tokenize(refseq_record, '\t', refseq_fields);
            size_t transcription_start = 0;
            if (refseq_fields[strand] == "+")
            {
                transcription_start =  refseq_fields[txStart].Int();
            }
            else if (refseq_fields[strand] == "-")
            {
                transcription_start = refseq_fields[txEnd].Int();
            }
            else
            {
                FatalErr("Unable to interpret strand field of " + 
                    refseq_record);
            }

            String chrom_name = refseq_fields[chrom].After("chr");
            
            if (chrom_name != "Y" && !chrom_name.Contains("_"))
            {
                int p = BinPosition(names, chrom_name);
                if (p < 0)
                {
                    FatalErr("Unable to locate " + chrom_name);
                }
                else
                {
                    int c = name_ids[p];
                    longlong wide_coverage = 0;
                    if (transcription_start < 1500 ||
                        transcription_start >= coverage[c].size() - 1500)
                    {
                        FatalErr(refseq_fields[name2] +
                                " wide window is off the reference");
                    }
                    for (size_t b = transcription_start - 1500;
                        b < transcription_start + 1500; b++)
                    {
                        wide_coverage += coverage[c][b];
                    }
                    longlong narrow_coverage = 0;
                    for (size_t b = transcription_start - 100;
                        b < transcription_start + 100; b++)
                    {
                        narrow_coverage += coverage[c][b];
                    }
                    double ratio = (narrow_coverage / 200.0) /
                        (wide_coverage / 3000.0);
                    bad_promoter_intervals.push_back(chrom_name + ":" + 
                        ToString(transcription_start - 100) + "-" +    
                        ToString(transcription_start + 100));
                    bad_promoter_ratios.push_back(ratio);
                    bad_promoter_genes.push_back(refseq_fields[name2]);
                }
            }
        }
    }
    
    // eliminate duplicate gene entries
    SortSync(bad_promoter_genes, bad_promoter_ratios, bad_promoter_intervals);
    vec<String> deduped_bad_promoter_genes;
    vec<double> deduped_bad_promoter_ratios;
    vec<String> deduped_bad_promoter_intervals;
    String current_gene;
    double current_gene_worst_ratio = -1;
    String current_interval;
    for (size_t b = 0; b < bad_promoter_genes.size(); b++)
    {
        if (current_gene != bad_promoter_genes[b])
        {
            current_gene = bad_promoter_genes[b];
            current_gene_worst_ratio = bad_promoter_ratios[b];
            current_interval = bad_promoter_intervals[b];
        }
        else if (bad_promoter_ratios[b] < current_gene_worst_ratio)
        {
            current_gene_worst_ratio = bad_promoter_ratios[b];
            current_interval = bad_promoter_intervals[b];
        }
        
        if (b == bad_promoter_genes.size() - 1 ||
            current_gene != bad_promoter_genes[b + 1])
        {
            deduped_bad_promoter_genes.push_back(current_gene);
            deduped_bad_promoter_ratios.push_back(current_gene_worst_ratio);
            deduped_bad_promoter_intervals.push_back(current_interval);
        }
    }
    
    // output the 1000 worst
    SortSync(deduped_bad_promoter_ratios, deduped_bad_promoter_intervals, 
        deduped_bad_promoter_genes);
    
    for (size_t b = 0; b < 1000 && b < deduped_bad_promoter_ratios.size(); b++)
    {
        std::cout << deduped_bad_promoter_intervals[b] << "\t#" << 
            deduped_bad_promoter_genes[b] << ", ratio=" << 
            deduped_bad_promoter_ratios[b] << "\n";
    }
    
    return 0;
}
