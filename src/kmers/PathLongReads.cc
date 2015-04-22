///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file PathLongReads.cc
 * \author tsharpe
 * \date Feb 28, 2012
 *
 * \brief
 */
#include "MainTools.h"
#include "feudal/VirtualMasterVec.h"
#include "kmers/LongReadPather.h"

#define LOGD(x) std::cout << Date() << ": " << x << std::endl;

int main( int argc, char** argv )
{
    // K must be odd:
    // we need to decide FWD or REV quickly,
    // and don't want to deal with palindromes
    unsigned const K = 10001;

    RunTime();
    BeginCommandArguments;
    CommandArgument_String_Doc(READS,"Fastb of reads to path.");
    CommandArgument_UnsignedInt_OrDefault_Doc(COVERAGE,30,"Estimated coverage.");
    CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS,0,"Max threads.");
    CommandArgument_String_OrDefault_Doc(OUTDIR,".","Output directory.");
    CommandArgument_Bool_OrDefault_Doc(WRITE_BINARY,False,"Binary versions of graph files.");
    CommandArgument_Bool_OrDefault_Doc(WRITE_FASTG,True,"FASTG output.");
#ifdef NDEBUG
    CommandArgument_UnsignedInt_OrDefault_Doc(LOGLEVEL,2,"Verbosity 0 to 5.");
#else
    CommandArgument_UnsignedInt_OrDefault_Doc(LOGLEVEL,5,"Verbosity 0 to 5.");
#endif
    EndCommandArguments;

    Long::gLogLevel = LOGLEVEL;

    Long::GraphFileNamer namer(OUTDIR,READS,K);
    String tmpDict(namer.getDictFilename() + ".tmp");
    VirtualMasterVec<bvec> reads(READS);
    LOGD("Building kmer dictionary.");
    Long::KmerDict<K> dict(4*reads.sizeSum(),COVERAGE);
    dict.addReads(reads.begin(),reads.end(),NUM_THREADS);
    std::cout << "There are " << dict.size() << " entries in the dictionary"
                << std::endl;

    LOGD("Building unipaths");
    Long::GraphBuilder<K>::EdgeVec edgeVec;
    Long::GraphBuilder<K> graphBuilder(dict,edgeVec);
    graphBuilder.buildUnipaths();

    LOGD("Joining unipaths.");
    graphBuilder.joinGraph();
    if ( WRITE_BINARY )
    {
        LOGD("Saving unipath graph.");
        graphBuilder.dumpGraph(namer);
    }
    if ( WRITE_FASTG )
    {
        LOGD("Writing FASTG.");
        graphBuilder.writeFASTG(namer);
    }
    graphBuilder.logStats();
    LOGD("Done.");
}
