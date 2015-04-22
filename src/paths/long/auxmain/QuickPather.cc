///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * QuickPather.cc
 *
 *  Created on: Nov 21, 2013
 *      Author: neilw
 */
#include "MainTools.h"
#include "FastaFileset.h"
#include "ParseSet.h"
#include "kmers/BigKPather.h"
#include <vector>

int main( int argc, char** argv )
{
    String empty;
    RunTime();
    BeginCommandArguments;
    CommandArgument_String_Doc(IN_FASTA,
        "looks for fasta files ParseStringSet format (i.e. multiple files)");
    CommandArgument_UnsignedInt_OrDefault_Doc(COVERAGE,30u,
        "Number of times each kmer occurs, on average, in fasta input.");
    CommandArgument_String_Doc(OUT_HEAD,
        "writes something");
    CommandArgument_UnsignedInt_OrDefault(NUM_THREADS,0u);
    CommandArgument_LongLong_OrDefault(MAX_MEMORY_GB,0);
    CommandArgument_UnsignedInt_OrDefault(K,60);
    EndCommandArguments;

    NUM_THREADS = configNumThreads(NUM_THREADS);
    SetMaxMemory(MAX_MEMORY_GB<<30);

    vecbasevector bases;

    cout << Date( ) << ": loading edges" << endl;
    vec<String> inputs;
    ParseStringSet( IN_FASTA, inputs );
    cout << Date( ) << ": found " << inputs.size( ) << " inputs" << endl;
    for ( int i = 0; i < inputs.isize( ); i++ ) {
        vecbasevector tmp_bases;
        FastFetchReads( tmp_bases, nullptr, inputs[i] );
        bases.resize( bases.size() + tmp_bases.size() );
        std::copy(tmp_bases.begin(), tmp_bases.end(),
                bases.end() - tmp_bases.size() );
    }

    cout << Date( ) << ": total edges = " << ToStringAddCommas( bases.size( ) )
         << endl;

    cout << Date() << ": pathing edges" << endl;

    HyperBasevector hb;
    buildBigKHBVFromReads(K,bases,COVERAGE,&hb);

    Ofstream( out, OUT_HEAD+".dot" );
    hb.PrintSummaryDOT0w( out, True, False, True );

    cout << Date() << ": Done." << std::endl;
}

