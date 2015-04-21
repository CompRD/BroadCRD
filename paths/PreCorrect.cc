///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file PreCorrect.cc
 * \author tsharpe
 * \date Nov 26, 2012
 *
 * \brief
 */

#include "MainTools.h"
#include "Basevector.h"
#include "bam/CachedBAMFile.h"
#include "paths/long/PreCorrectAlt1.h"
#include "simulation/BamReadGenerator.h"
#include "simulation/ReadErrorAnalyzer.h"
#include <algorithm>
#include <vector>

namespace
{
    void buildRefReads( String const& refFastb, RefLocusVec const& refLocs,
                            vecbvec const& reads, vecbvec* pVBV )
    {
        pVBV->reserve(refLocs.size());
        vecbvec ref(refFastb);
        bvec tmp;
        typedef RefLocusVec::const_iterator Itr;
        vecbvec::const_iterator rdItr(reads.begin());
        for ( Itr itr(refLocs.begin()), end(refLocs.end()); itr != end; ++itr )
        {
            RefLocus const& loc = *itr;
            bvec const& refBV = ref[loc.mRefID];
            bvec::const_iterator beg(refBV.begin(loc.mOffset));
            size_t readLen = rdItr->size();
            ++rdItr;
            int wrap = loc.mOffset + readLen - refBV.size();
            if ( wrap <= 0 )
                tmp.assign(beg,beg+readLen);
            else
                tmp.assign(beg,refBV.end()).append(refBV.begin(),refBV.begin(wrap));
            if ( loc.mRC )
                tmp.ReverseComplement();
            pVBV->push_back(tmp);
        }
    }
}

int main( int argc, char** argv )
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_UnsignedInt_OrDefault_Doc(COVERAGE,50u,
            "Estimated coverage.");
    CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS,0u,
            "Number of parallel threads.");
    CommandArgument_Int_OrDefault_Doc(VERBOSITY,0,
            "Document each correction.");
    EndCommandArguments;

    NUM_THREADS = configNumThreads(NUM_THREADS);

    std::vector<CachedBAMFile> bams;
    bams.push_back(CachedBAMFile("A1PJD","C1-508_2012-09-24_2012-10-04",1,"Solexa-121745.Tag_IlluminaHTKit2"));
    bams.push_back(CachedBAMFile("A1PJD","C1-508_2012-09-24_2012-10-04",1,"Solexa-121745.Tag_IlluminaHTKit3"));
    bams.push_back(CachedBAMFile("A1PJD","C1-508_2012-09-24_2012-10-04",1,"Solexa-121745.Tag_IlluminaHTKit4"));
    bams.push_back(CachedBAMFile("A1PJD","C1-508_2012-09-24_2012-10-04",1,"Solexa-121745.Tag_IlluminaHTKit5"));

    std::cout << Date() << ": Reading BAMs." << std::endl;
    BamReadGenerator brg(bams);
    vecbvec reads(brg.getReads());
    if ( reads.empty() )
        FatalErr("There are no reads to process.");

    vecbvec corrReads(reads);
    precorrectAlt1(&corrReads,COVERAGE,VERBOSITY,NUM_THREADS);

    std::cout << Date() << " Evaluating corrections." << std::endl;
    vecbvec refReads;
    buildRefReads("/wga/dev/references/Rhodobacter_sphaeroides/genome.fastb",
                    brg.getReadLocs(),reads,&refReads);

    ReadErrorAnalyzer rea(refReads, reads, corrReads);
    rea.SetVerbosity(VERBOSITY);
    rea.Analyze();
    std::cout << Date() << ": Done." << std::endl;
}
