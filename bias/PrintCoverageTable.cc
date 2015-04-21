///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Print coverage in a specified interval.

#include "MainTools.h"
#include "bias/CoverageIterator.h"
#include "bias/GenomeReference.h"

int main(int argc, char** argv)
{
    RunTime();
    
    BeginCommandArguments;
    CommandArgument_String_Doc(COV, "coverage BAM/SCI file");
    CommandArgument_String_Doc(REF, "genome reference file");
    CommandArgument_String_Doc(INTERVAL, "genome location to print");
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

    String contig = INTERVAL.Before(":");
    String range = INTERVAL.After(":");
    size_t start = range.Before("-").Int();
    size_t stop = range.After("-").Int();
    bool interval_printed = false;
    const vec<String>& names = gref.getNames();
    const vecbitvector& refamb = gref.getAmbiguities();
    const vecbasevector& genome = gref.getBases();
    
    while (coviter->hasNext() && !interval_printed)
    {
        BC::Coverage c = coviter->next();
        if (names[c.chr] == contig)
        {
            if (c.loc >= start && c.loc < stop)
            {
                std::cout << contig << ":" << c.loc << "\t";
                if (!refamb[c.chr][c.loc])
                { 
                    std::cout << Base::val2Char(genome[c.chr][c.loc]) << "\t";
                }
                else
                {
                    std::cout << "N\t";
                }
                std::cout << (c.fwd + c.rev) << "\n";
            }
            else if (c.loc >= stop)
            {
                interval_printed = true;
            }
        }
    }

    return 0;
}