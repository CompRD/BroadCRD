/*
 * KmerizeAndVote.cc
 *
 *  Created on: Dec 19, 2013
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

#include "Histogram.h"

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

    CommandArgument_String_Doc(IN_HEAD, "Looks for IN_HEAD.fasta for the consensus");
    CommandArgument_String_Doc(REF, "Looks for REF.fasta for reference");
    CommandArgument_String_Doc(BANK_HEAD, "bank of reads in BANK_HEAD.{fastb,qualb}");
    CommandArgument_String_Doc(RANK, "Looks for RANK for list of rank and indices ");
    CommandArgument_UnsignedInt_Doc(NREADS, "maximum number of reads to participate");
//    CommandArgument_UnsignedInt_Doc(K, "K value");
    CommandArgument_String_OrDefault_Doc(MOTIF, "", "generate in depth analysis for a particular motif");
    CommandArgument_UnsignedInt_Doc(THRESHOLD, "quality-score threshold");
    CommandArgument_UnsignedInt_OrDefault_Doc(POS_DEV, 1, "position uncertainty");
    CommandArgument_UnsignedInt_OrDefault_Doc(FLANK, 1, "number of flanking base of a certain motif");

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
    vecqualvector quals;
    if(IsRegularFile(RANK)){
        std::cout << Date() << ": loading read bank" << BANK_HEAD << ".fastb" << std::endl;
        vecbasevector read_bank(BANK_HEAD+".fastb");
        std::cout << Date() << ": loading read bank" << BANK_HEAD << ".qualb" << std::endl;
        vecqualvector qual_bank(BANK_HEAD+".qualb");
        ForceAssert(read_bank.size()==qual_bank.size());
        std::cout << Date() << ": number of reads in bank: " << read_bank.size() << std::endl;
        std::cout << Date() << ": loading rank file " << RANK << std::endl;
        ifstream ifs(RANK);
        int last_score=-1;
        int i,j;
        for(ifs>>i>>j ; reads.size()<NREADS && reads.size() < read_bank.size() && ifs.good() ; ifs>>i>>j){
            ForceAssert( i >= last_score);
            last_score=i;
            reads.push_back( read_bank[j] );
            quals.push_back( qual_bank[j] );
        }
        ifs.close();
        std::cout << Date() << ": done" << std::endl;
    }
    else{
        FatalErr("Cannot open rank file.");
    }

    std::cout << Date() << ": Number of reads " << reads.size() << std::endl;
    for(auto& read:reads){
        read.SetToSubOf(read,can_begin,-1);
        longest = max(longest,read.size());
    }

    vec<alignment> read_alignments(reads.size());
    vec<vec<size_t>> map12(reads.size());
    std::cout << Date() << ": aligning reads to reference " << std::endl;
    #pragma omp parallel
    {
        pacbio_bridge_tools::SmithWatAffineResource_t resource(longest+1,longest+1);
        #pragma omp for
        for(size_t rr=0;rr<reads.size();++rr){
            pacbio_bridge_tools::SmithWatAffine_loc(reads[rr],reference[0],read_alignments[rr],resource, false,false,1,1,1);
            pacbio_bridge_tools::pos1_to_pos2(map12[rr],read_alignments[rr]);
            ForceAssert(map12[rr].size()==reads[rr].size());
    //        PrintVisualAlignment(True,std::cout,reads[rr],reference[0],read_alignments[rr]);
        }
    }

    if(MOTIF!=""){
        basevector bMotif(MOTIF);
        std::cout << Date() << ": motif of length " << bMotif.size() << " : " << bMotif.ToString() << std::endl;
        pacbio_bridge_tools::kp_analysis_t kp_analysis(bMotif.size());
        kp_analysis.processReads(reads,quals,map12,THRESHOLD,5);
        std::cout << Date() << ": Report:" << std::endl;
        kp_analysis.printKmerDistribution( pacbio_bridge_tools::KmerVal(bMotif,bMotif.size()) );
    }

    pacbio_bridge_tools::SmithWatAffineResource_t resource(longest+1,longest+1);
    for(size_t cc=0 ; cc < candidates.size() ; ++cc){
        const auto& candidate = candidates[cc];
        alignment a;
        pacbio_bridge_tools::SmithWatAffine_loc(candidate,reference[0],a,resource, false,false,1,1,1);
        std::cout << "alignment of consensus " << cc << std::endl;
        PrintVisualAlignment(True,std::cout,candidate,reference[0],a);
        pacbio_bridge_tools::alternatives_t alt(candidate,reference[0],a,FLANK);
        alt.pileup(reads,quals,map12,THRESHOLD,FLANK,POS_DEV);
        std::cout << "pileup results" << std::endl;
        std::cout << alt << std::endl;
    }

    std::cout << Date() << ": Exiting." << std::endl;
}
