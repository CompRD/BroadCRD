///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/* This program measures the correlation in coverage, quality, and error rates
   across two BAM files. */

// Original author: Michael G. Ross <mgross@broadinstitute.org>

#include "MainTools.h"
#include "bias/CoverageIterator.h"
#include "bias/GenomeReference.h"

int main(int argc, char** argv)
{
    RunTime();
    
    BeginCommandArguments;
    CommandArgument_String_Doc(COV1, "first coverage BAM/SCI file");
    CommandArgument_String_Doc(COV2, "second coverage BAM/SCI file");
    CommandArgument_String_Doc(REF, "genome reference file");
    CommandArgument_String_OrDefault_Doc(REF_CACHE_DIR, "",
        "directory of cached genome reference information");
    EndCommandArguments;
    
    REF = RealPath(REF);
    if (!REF_CACHE_DIR.empty())
    {
        REF_CACHE_DIR = REF_CACHE_DIR + "/" + REF;
    }
    
    BC::GenomeReference gref(REF, REF_CACHE_DIR);
    
    BC::AbstractCoverageIterator* coviter1 = load_coverage_iter(COV1, gref);
    BC::AbstractCoverageIterator* coviter2 = load_coverage_iter(COV2, gref);
    
    longlong cov1_total = 0;
    longlong cov2_total = 0;
    longlong cov1sq_total = 0;
    longlong cov2sq_total = 0;
    longlong cov1cov2_total = 0;
    longlong gen_size = 0;
    
    while (coviter1->hasNext() && coviter2->hasNext())
    {
        BC::Coverage c1 = coviter1->next();
        BC::Coverage c2 = coviter2->next();
        
        longlong cov1 = c1.fwd + c1.rev;
        longlong cov2 = c2.fwd + c2.rev;
        
        cov1_total += cov1;
        cov2_total += cov2;
        cov1sq_total += cov1 * cov1;
        cov2sq_total += cov2 * cov2;
        cov1cov2_total += cov1 * cov2;
        
        gen_size++;
    }
    
    if (coviter1->hasNext() || coviter2->hasNext())
    {
        std::cerr << "somehow the iterators are not properly synchronized"
            << std::endl;
        return 1;
    }
    
    double cov1_mean = static_cast<double>(cov1_total) / gen_size;
    double cov2_mean = static_cast<double>(cov2_total) / gen_size;
    double cov1sq_mean = static_cast<double>(cov1sq_total) / gen_size;
    double cov2sq_mean = static_cast<double>(cov2sq_total) / gen_size;
    double cov1cov2_mean = static_cast<double>(cov1cov2_total) / gen_size;

    PRINT5(cov1_mean, cov2_mean, cov1sq_mean, cov2sq_mean, cov1cov2_mean);

    double cov1_var = cov1sq_mean - pow(cov1_mean, 2);
    double cov2_var = cov2sq_mean - pow(cov2_mean, 2);
    double cov1cov2_covar = cov1cov2_mean - cov1_mean * cov2_mean;
    
    PRINT3(cov1_var, cov2_var, cov1cov2_covar);
    
    double correlation = cov1cov2_covar / (sqrt(cov1_var) * sqrt(cov2_var));
    
    std::cout << "coverage1 variance = " << cov1_var << std::endl;
    std::cout << "coverage2 variance = " << cov2_var << std::endl;
    std::cout << "coverage correlation = " << correlation << std::endl;
    
    delete coviter1;
    delete coviter2;
    
    return 0;
}
