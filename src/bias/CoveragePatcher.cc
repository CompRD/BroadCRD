///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// This program reports statistics indicating how well coverage from one
// technology can patch holes in the coverage from another technology.

// Original author: Michael G. Ross <mgross@broadinstitute.org>

#include "MainTools.h"
#include "Vec.h"
#include "bias/CoverageIterator.h"
#include "bias/GenomeReference.h"

vec<longlong> accumulate_counts(vec<longlong> counts)
{
    vec<longlong> accumulated(counts.size(), 0);
    accumulated[0] = counts[0];
    for (size_t i = 1; i < counts.size(); i++)
    {
        accumulated[i] = accumulated[i - 1] + counts[i];
    }
    return accumulated;
}

int main(int argc, char** argv)
{
    RunTime();
    
    BeginCommandArguments;
    CommandArgument_String_Doc(ORIGCOV, "original coverage BAM/SCI files");
    CommandArgument_String_Doc(PATCHCOV, "patching coverage BAM/SCI files");
    CommandArgument_String_OrDefault_Doc(REF, "", "genome reference file");
    CommandArgument_String_OrDefault_Doc(REF_CACHE_DIR, "",
        "directory of cached genome reference information");
    EndCommandArguments;

    vec<String> orig_files = BC::locate_cov_files(ORIGCOV);
    vec<String> patch_files = BC::locate_cov_files(PATCHCOV);

    std::cout << "original coverage files = ";
    CompactPrint(std::cout, orig_files, ", ");
    std::cout << "\n";
    std::cout << "patch coverage files = ";
    CompactPrint(std::cout, patch_files, ", ");
    std::cout << "\n";

    if (REF.empty())
    {
        String ref1 = BC::locate_reference(orig_files);
        String ref2 = BC::locate_reference(patch_files);
        
        if (ref1.empty())
        {
            FatalErr("Unable to infer original reference.");
        }
        else if (ref2.empty())
        {
            FatalErr("Unable to infer patching reference.");
        }
        
        if (ref1 == ref2)
        {
            REF = ref1;
        }
        else
        {
            FatalErr("Original and patch files are using different references: "
                + ref1 + " and " + ref2 + ".");
        }
    }
    
    REF = RealPath(REF);

    BC::GenomeReference gref = BC::GenomeReference(REF, REF_CACHE_DIR + "/" +
        REF);

    BC::MultiCoverageIterator orig_coviter(orig_files, gref);
    BC::MultiCoverageIterator patch_coviter(patch_files, gref);

    longlong mtotal_orig_cov = 0;
    longlong mtotal_patch_cov = 0;
    longlong loc_count = 0;
    while (orig_coviter.hasNext())
    {
        BC::Coverage co = orig_coviter.next();
        BC::Coverage cp = patch_coviter.next();
        if (!gref.getAmbiguities()[co.chr][co.loc])
        {
            mtotal_orig_cov += co.fwd + co.rev;
            mtotal_patch_cov += cp.fwd + cp.rev;            
            loc_count++;
        }
    }

    double mean_orig_cov = static_cast<double>(mtotal_orig_cov) / loc_count;
    double mean_patch_cov = static_cast<double>(mtotal_patch_cov) / loc_count;

    std::cout << "mean_orig_cov = " << mean_orig_cov << "\n";
    std::cout << "mean_patch_cov = " << mean_patch_cov << "\n";
    
    const int NUM_BINS = static_cast<int>(mean_orig_cov);
            
    vec<longlong> total_patch_cov(NUM_BINS, 0);
    vec<longlong> total_orig_cov(NUM_BINS, 0);
    vec<longlong> total_patch_err(NUM_BINS, 0);
    vec<longlong> total_orig_err(NUM_BINS, 0);
    vec<longlong> patch_counts(NUM_BINS, 0);
    
    orig_coviter.rewind();
    patch_coviter.rewind();
    while (orig_coviter.hasNext())
    {
        BC::Coverage co = orig_coviter.next();
        BC::Coverage cp = patch_coviter.next();
        
        if (!gref.getAmbiguities()[co.chr][co.loc])
        {
            longlong orig_cov = co.fwd + co.rev;
            longlong patch_cov = cp.fwd + cp.rev;
            longlong orig_err = co.fwd_mismatches + co.rev_mismatches
                + co.fwd_insertions + co.rev_insertions
                + co.fwd_deletions + co.rev_deletions;
            longlong patch_err = cp.fwd_mismatches + cp.rev_mismatches
                + cp.fwd_insertions + cp.rev_insertions
                + cp.fwd_deletions + cp.rev_deletions;
    
            double orig_rcov = orig_cov / mean_orig_cov;
            
            int slot = static_cast<int>(orig_rcov * NUM_BINS);
            if (slot < NUM_BINS)
            {
                total_orig_cov[slot] += orig_cov;
                total_patch_cov[slot] += patch_cov;
                total_orig_err[slot] += orig_err;
                total_patch_err[slot] += patch_err;
                patch_counts[slot]++;
            }
        }
    }
    
    vec<longlong> accum_orig_cov = accumulate_counts(total_orig_cov);
    vec<longlong> accum_patch_cov = accumulate_counts(total_patch_cov);
    vec<longlong> accum_orig_err = accumulate_counts(total_orig_err);
    vec<longlong> accum_patch_err = accumulate_counts(total_patch_err);
    vec<longlong> accum_patch_counts = accumulate_counts(patch_counts);
        
    vec<double> mean_patch_rcov(NUM_BINS, 0);
    vec<double> mean_orig_err(NUM_BINS, 0);
    vec<double> mean_patch_err(NUM_BINS, 0);
    vec<double> accum_mean_patch_rcov(NUM_BINS, 0);
    vec<double> accum_mean_orig_err(NUM_BINS, 0);
    vec<double> accum_mean_patch_err(NUM_BINS, 0);

    for (size_t i = 0; i < total_orig_cov.size(); i++)
    {
        mean_patch_rcov[i] = static_cast<double>(total_patch_cov[i]) / 
            (patch_counts[i] * mean_patch_cov);
        mean_orig_err[i] = static_cast<double>(total_orig_err[i]) / 
            total_orig_cov[i];
        mean_patch_err[i] = static_cast<double>(total_patch_err[i]) /
            total_patch_cov[i];
        accum_mean_patch_rcov[i] = static_cast<double>(accum_patch_cov[i])/
            (accum_patch_counts[i] * mean_patch_cov);
        accum_mean_orig_err[i] = static_cast<double>(accum_orig_err[i]) /
            accum_orig_cov[i];
        accum_mean_patch_err[i] = static_cast<double>(accum_patch_err[i]) /
            accum_patch_cov[i];
    }

    std::cout << "\n";
    std::cout << "original rel coverage\tmean patch rel coverage\t"
        "mean orig erors\tmean patch errors\tpatch counts\n";
    for (size_t i = 0; i < total_orig_cov.size(); i++)
    {
        std::cout << ((i + 1) / static_cast<double>(NUM_BINS)) << "\t" 
            << mean_patch_rcov[i] << "\t" << mean_orig_err[i] << "\t"
            << mean_patch_err[i] << "\t" << patch_counts[i] << "\n";
    }

    std::cout << "\n\n";
    std::cout << "original rel coverage\taccum mean patch rel coverage\t"
        "accum mean orig errors\taccum mean patch errors\t"
        "accum patch counts\n";
    for (size_t i = 0; i < total_orig_cov.size(); i++)
    {
        std::cout << ((i + 1) / static_cast<double>(NUM_BINS)) << "\t"
            << accum_mean_patch_rcov[i] << "\t" << accum_mean_orig_err[i]
            << "\t" << accum_mean_patch_err[i] << "\t"
            << accum_patch_counts[i] << "\n";
    }

    return 0;
}
