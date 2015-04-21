/*
 * PileUpAndVote.cc
 *
 *  Created on: Dec 18, 2013
 *      Author: blau
 */

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
#include <limits>
#include <fstream>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <unordered_map>
#include <functional>
#include <iterator>
#include <omp.h>

#include "Basevector.h"
#include "FetchReads.h"
#include "PrintAlignment.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/LongReadPatchOptimizer.h"
#include "pacbio/BridgeTools.h"


#include "MainTools.h"
int main(int argc, char** argv) {
    RunTime();
    BeginCommandArguments;

    CommandArgument_String_Doc(IN_HEAD, "Looks for IN_HEAD.fasta for the reads");
    CommandArgument_String_Doc(REF, "Looks for REF.fasta for reference");
    CommandArgument_String_Doc(READ_HEAD, "Looks for READ_HEAD.fasta for reference");

    EndCommandArguments;

    vecbasevector reference;
    FetchReads(reference, 0, REF+".fasta");
    const size_t ref_begin = 10;
    reference[0].SetToSubOf(reference[0],ref_begin,-1);
    unsigned int longest=reference[0].size();

    vecbasevector candidates;
    FetchReads(candidates, 0, IN_HEAD+".fasta");
    const size_t can_begin = 37;
    for(auto& candidate : candidates){
        candidate.SetToSubOf(candidate,can_begin,-1);
        longest = max(longest,candidate.size());
    }
    std::cout << Date() << ": Number of candidates " << candidates.size() << std::endl;

    vecbasevector reads;
    FetchReads(reads, 0, READ_HEAD+".fasta");
    for(auto& read:reads){
        read.SetToSubOf(read,can_begin,-1);
        longest = max(longest,read.size());
    }
    vec<alignment> read_alignments(reads.size());
    std::cout << Date() << ": Number of reads " << reads.size() << std::endl;



    #pragma omp parallel
    {
        pacbio_bridge_tools::SmithWatAffineResource_t resource(longest+1,longest+1);
        #pragma omp for
        for(size_t rr=0;rr<reads.size();++rr){
    //        std::cout << "alignment of read " << rr << std::endl;
            pacbio_bridge_tools::SmithWatAffine_loc(reads[rr],reference[0],read_alignments[rr],resource, false,false,1,1,1);
    //        PrintVisualAlignment(True,std::cout,reads[rr],reference[0],read_alignments[rr]);
        }
    }


    pacbio_bridge_tools::SmithWatAffineResource_t resource(longest+1,longest+1);
    for(size_t cc=0 ; cc < candidates.size() ; ++cc){
        const auto& candidate = candidates[cc];
        alignment a;
        pacbio_bridge_tools::SmithWatAffine_loc(candidate,reference[0],a,resource, false,false,1,1,1);
        std::cout << "alignment of consensus " << cc << std::endl;
        PrintVisualAlignment(True,std::cout,candidate,reference[0],a);
        pacbio_bridge_tools::alternatives_t alt(candidate,reference[0],a);
        alt.pileup(reads,read_alignments);
        std::cout << "pileup results" << std::endl;
        std::cout << alt << std::endl;
    }

    std::cout << Date() << ": Done." << std::endl;
}
