///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file TestMultipleAligner.cc
 * \author tsharpe
 * \date May 11, 2012
 *
 * \brief
 */
#include "MainTools.h"
#include "paths/long/ultra/MultipleAligner.h"
#include <iostream>
#include <fstream>

int main( int argc, char** argv )
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_Double_OrDefault_Doc(ERR_DEL, 0.03, "deletion rate");
    CommandArgument_Double_OrDefault_Doc(ERR_INS, 0.002, "insertion rate");
    CommandArgument_Double_OrDefault_Doc(ERR_SUB, 0.008, "substitution rate");
    CommandArgument_String_Doc(CONSENSUS,"consensus base calls");
    CommandArgument_String_Doc(THREADS_FASTA,"fasta file for thread calls");
    EndCommandArguments;

    vecbvec reads;
    std::fstream fasta(THREADS_FASTA.c_str());
    std::string line;
    while ( getline(fasta,line) )
        if ( line.size() && line[0] != '>' )
            reads.push_back(bvec(line));

    Scorer scorer( ERR_SUB, ERR_DEL, ERR_INS );
    bvec consensus(CONSENSUS);
    MultipleAligner ma(scorer,consensus);
    ma.addReads(reads.begin(),reads.end());
    ma.printMultipleAlignment(std::cout,reads.begin(),reads.end());
}
