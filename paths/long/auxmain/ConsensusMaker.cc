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
#include "pairwise_aligners/SmithWatAffine.h"
#include "PrintAlignment.h"
#include "lookup/LookAlign.h"


int main(int argc, char* argv[])
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_Int_OrDefault(K, 21);
    CommandArgument_Int_OrDefault_Doc(NUM_THREADS, 0,
	    "number of threads to use in kmerization");
    CommandArgument_String_Doc(IN_HEAD, "head for reads to process");
    CommandArgument_String_Doc(OUT_HEAD, "output head for files");
    CommandArgument_UnsignedInt_OrDefault_Doc(MAX_FREQ, 1000, "maximum kmer frequency allowed" );
    CommandArgument_String_OrDefault_Doc(TMP_HEAD, "./tempxxx", "head for large, temporary files" );
    EndCommandArguments;

    NUM_THREADS = configNumThreads(NUM_THREADS);
    omp_set_num_threads(NUM_THREADS);

    // read input files
    cout << Date() << ": reading in reads" << endl;
    BaseVecVec readvv(IN_HEAD + ".fastb");

    // create alignments
    vec<look_align> aligns;
    vec< vec<int> > aligns_index;
    LoadLookAligns( IN_HEAD+".aligns", aligns, aligns_index, readvv.size( ) );

    for ( size_t readid = 0; readid < readvv.size(); ++readid ) {
	auto const& founder = readvv[readid];

	for ( auto const& friendalignid : aligns_index[readid] ) {
	    auto const& align = aligns[friendalignid];
	    PRINT3( readid, align.query_id, align.target_id );

	    auto const& buddy = readvv[align.target_id];

	    alignment a;
	    SmithWatAffine(buddy,founder,a,true,true);

	    PrintVisualAlignment(true,cout,buddy,founder,a);
	}
    }


#if 0
    typedef Kmer60 Kmer_t;

    cout << Date() << ": about to instantiate the FriendAlignFinderNaif" << endl;
    FriendAlignFinderNaif<Kmer_t> ff( TMP_HEAD, K, readvv, MAX_FREQ, NUM_THREADS );

    for ( size_t readid = 0; readid < readvv.size(); ++readid ) {
	auto const& founder = readvv[readid];

	Friends f;
	ff.getAligns(readid, &f);

	cout << "readid=" << readid << ", friends=" << f.size() << endl;

	for ( size_t friendid = 0; friendid < f.size(); ++friendid ) {
	    size_t friend_readid = f[friendid].readId();
	    auto const& buddy = readvv[friend_readid];
	    alignment a;
//	    SmithWatAffineSuper(buddy,founder,a,61,false,false,0,3,12,1,SmithWatAffineParallel);
	    SmithWatAffine(buddy,founder,a,false,false);

	    PrintVisualAlignment(true,cout,buddy,founder,a);
	}
    }

#endif

    cout << Date() << ": about to exit" << endl;

    return 0;

}

