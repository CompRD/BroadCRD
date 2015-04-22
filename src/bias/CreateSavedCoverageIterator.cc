///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/* This program converts multiple BAM/SAM/SCI files into a single
   SavedCoverageIterator (SCI) file, which can then be reprocessed more quickly
   by BadCoverage or other coverage analysis tools. */

// Original author: Michael G. Ross <mgross@broadinstitute.org>

#include "MainTools.h"
#include "bias/CoverageIterator.h"
#include "bias/GenomeReference.h"

int main(int argc, char** argv)
{
    RunTime();
    
    BeginCommandArguments;
    CommandArgument_String_Doc(COV, "coverage files");
    CommandArgument_String_Doc(REF, "genome reference file");
    CommandArgument_String_Doc(SCI, "output SCI file");
    CommandArgument_String_OrDefault_Doc(REF_CACHE_DIR, "",
        "directory of cached genome reference information");
    EndCommandArguments;

    REF = RealPath(REF);
    if (!REF_CACHE_DIR.empty())
    {
        REF_CACHE_DIR = REF_CACHE_DIR + "/" + REF;
    }

    BC::GenomeReference gref(REF, REF_CACHE_DIR);
    
    vec<String> cov_files = BC::locate_cov_files(COV);
    
    std::cout << "input files = ";
    CompactPrint(std::cout, cov_files, ", ");
    std::cout << "\n";
    
    BC::MultiCoverageIterator coviter(cov_files, gref);
    BC::SavedCoverageIterator sciter(SCI, REF, coviter.getLibraries(),
        coviter.getProductionUnits(), coviter.getSamples());
        
    while (coviter.hasNext())
    {
        sciter.addCoverage(coviter.next(), coviter.numReads(),
            coviter.numReadsUsed());
    }
    
    return 0;
}
