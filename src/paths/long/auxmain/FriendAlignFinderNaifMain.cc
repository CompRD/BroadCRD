///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Nov 6, 2012
//

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "paths/long/FriendAlignFinderNaif.h"
#include "feudal/FilesOutputIterator.h"
#include "feudal/BinaryStream.h"
#include "ParallelVecUtilities.h"
#include "Qualvector.h"

namespace {

}

int main(int argc, char* argv[])
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_Int_OrDefault(K, 21);
    CommandArgument_Int_OrDefault_Doc(NUM_THREADS, 0,
	    "number of threads to use in kmerization");
    CommandArgument_String_Doc(FASTB, "Reads to process");
    CommandArgument_String_OrDefault_Doc(QUALB, "", "Quality scores to process");
    CommandArgument_String_Doc(OUT_HEAD, "output head for files");
    CommandArgument_UnsignedInt_OrDefault_Doc(MAX_FREQ, 1000, "maximum kmer frequency allowed" );
    CommandArgument_String_OrDefault_Doc(TMP_DIR, ".", "directory for large, temporary files" );
    EndCommandArguments;

    NUM_THREADS = configNumThreads(NUM_THREADS);
    omp_set_num_threads(NUM_THREADS);

    // read input files
    cout << Date() << ": reading in reads" << endl;
    BaseVecVec readvv(FASTB);

    // trim reads for Q20
    if ( QUALB != "" ) {
	QualVecVec qvv(QUALB);
	ForceAssertEq(readvv.size(), qvv.size() );
#pragma omp parallel for
	for ( size_t i = 0; i < readvv.size(); ++i ) {
	    ForceAssertEq(readvv[i].size(), qvv[i].size());
	    size_t trim = 0;
	    while ( trim < qvv[i].size() ) {
		if ( qvv[i][trim] < 20 ) break;
		trim++;
	    }
	    readvv[i].resize(trim);
	    qvv[i].resize(trim);
	}

    }




    // create alignments

    typedef Kmer60 Kmer_t;

    cout << Date() << ": about to instantiate the FriendAlignFinderNaif" << endl;
    FriendAlignFinderNaif<Kmer_t> ff( TMP_DIR, K, readvv, MAX_FREQ, NUM_THREADS );


    cout << Date() << ": about to exit" << endl;

    return 0;

}

