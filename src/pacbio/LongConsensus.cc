/*
 * LongConsensus.cc
 *
 *  Created on: Dec 13, 2013
 *      Author: blau
 */

#include <iostream>
#include "Basevector.h"
#include "paths/LongReadPatchOptimizer.h"



#include "MainTools.h"
int main(int argc, char** argv) {
    RunTime();
    BeginCommandArguments;

    CommandArgument_String_Doc(IN_HEAD, "Looks for IN_HEAD.fastb for the reads");
    CommandArgument_String_Doc(GUESS_HEAD, "Initial guess in GUESS_HEAD.fastb ");

    CommandArgument_Int_OrDefault_Doc(VERBOSITY,0,"verbosity");
    CommandArgument_Int_OrDefault_Doc(PAD,1,"padding");

    vecbasevector reads(IN_HEAD+".fastb");
    const String OUT_HEAD=IN_HEAD+".consensus";

    vecbasevector consensus(GUESS_HEAD+".fastb");
    if(consensus.size()>1){
        std::cout<<"WARNING: only the first sequence in " << GUESS_HEAD << ".fastb is used." << std::endl;
    }

    if(VERBOSITY){
        std::cout << ">Input sequences:\n";
        for(const auto& read: reads){
            std::cout << ">nBases=" << read.size() << "\n";
            std::cout << read.ToString() << std::endl;
        }
        std::cout << std::endl;

        std::cout << ">Initial Consensus\n";
        std::cout << consensus.front().ToString() << std::endl;
        std::cout << std::endl;
    }

    std::cout<<Date()<<"calling consensus_compute with " << reads.size() << " reads"<<std::endl;

    consensus_compute(reads,&consensus.front(),PAD,VERBOSITY);

    if(VERBOSITY){
        std::cout << std::endl;
        std::cout << ">Final Consensus\n";
        std::cout << consensus.front().ToString() << std::endl;
        std::cout << std::endl;
    }
    consensus.WriteAll(OUT_HEAD+".fastb");
    std::cout << Date() << ": Done." << std::endl;
}
