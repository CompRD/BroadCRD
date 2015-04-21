///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "Basevector.h"
#include "Bitvector.h"
#include "FeudalMimic.h"
#include "FetchReads.h"
#include "system/file/File.h"
#include "paths/long/MakeKmerStuff.h"
#include "system/System.h"
#include "kmers/naif_kmer/Kmers.h"
#include "kmers/naif_kmer/KernelKmerStorer.h"
#include "kmers/naif_kmer/KernelKmerMasker.h"

#undef COMPRD		// define COMPRD for additional options

static const int gK = 100;


void naif_build_repeat_mask( const String& REF, const String& dest, const size_t NUM_THREADS )
{
    cout << Date() << ": reading the reference" << endl;
    vecbasevector X;
    if ( REF.EndsWith(".fastb") ) X.ReadAll( REF );	// secret shortcut -- we can have a .fastb file
    else FetchReads( X, 0, REF );

    cout << Date() << ": kmerizing..." << endl;

    typedef Kmer124 Kmer_t;
    typedef KmerGVLoc<Kmer_t> KmerRec_t;

    vecbitvector bitvec;

    // calculate the kmer frequency db, thresholding at min_freq
    Validator valid( 2, 0 );	// only mark kmers >= 2
    KernelKmerMasker<KmerRec_t> storer( X, gK, &bitvec, &valid);
    naif_kmerize( &storer, NUM_THREADS, 1 );

    bitvec.WriteAll(dest);
}

#ifdef COMPRD
void naive_build_repeat_mask( const String& REF, const String& dest)
{
    vecbasevector X;
    FetchReads( X, 0, REF );
    vecbitvector multi;
    Mimic(X, multi);
    vecbasevector X_rc(X);
    for (int g = 0; g < (int) X.size(); g++)
	X_rc[g].ReverseComplement();
    X.Append(X_rc);

    vec<triple<kmer<gK>, int, int> > kmers_plus;
    cout << Date() << ": making kmer lookup" << endl;
    MakeKmerLookup1(X, kmers_plus);


    cout << Date() << ": looking for multis" << endl;
    for (int64_t i = 0; i < kmers_plus.jsize(); i++) {
	int64_t j;
	for (j = i + 1; j < kmers_plus.jsize(); j++)
	    if (kmers_plus[j].first != kmers_plus[i].first)
		break;
	if (j - i >= 2) {
	    for (int64_t k = i; k < j; k++) {
		int g = kmers_plus[k].second, p = kmers_plus[k].third;
		if (g < (int) X.size() / 2)
		    multi[g].Set(p, True);
		else {
		    g -= (int) X.size() / 2;
		    p = X[g].isize() - p - gK;
		    multi[g].Set(p, True);
		}
	    }
	}
	i = j - 1;
    }
    cout << Date() << ": writing " << dest << endl;

    multi.WriteAll(dest);
}
#endif


void build_index_if_needed( const String& ref )
{
    const String index=ref+".fai";
    const String samtools = "samtools";

    cout << "Checking fasta index " << index  << "... ";
    if ( File(index).exists() )  {
	cout << "found." << endl;
    } else {
	cout << "not found." << endl;
	cout << Date() << ": The fasta index does not exist, so we're going to build it..." << endl;
	if ( System("which "+samtools) != 0 ) FatalErr("Discovar requires that samtools be installed and in your path");
	SystemSucceed( samtools + " faidx " + ref );
	cout << Date() << ": done building index " + index << endl;
    }
}


void build_repeat_mask_if_needed( const String& REF, const size_t NUM_THREADS)
{
    String dest;
    String newExt = ".m"+ToString(gK);

    dest = REF + newExt;

    cout << "Checking reference mask " << dest << "... ";
    if ( File(dest).exists() )  {
	cout << "found." << endl;
    } else {
	cout << "not found." << endl;
	cout << Date() << ": The reference mask does not exist, so we're going to build it.  This will take a while..." << endl;
	naif_build_repeat_mask(REF, dest, NUM_THREADS);
//	naive_build_repeat_mask(REF, dest );
	cout << Date() << ": Done building reference mask " + dest << endl;
    }
}



int main(int argc, char* argv[] )
{
    RunTime( );

    BeginCommandArguments;
    CommandArgument_String_Doc( REF, "reference file in FASTA format");
    CommandArgument_UnsignedInt_OrDefault( NUM_THREADS, 0 );
    EndCommandArguments;

    NUM_THREADS = configNumThreads(NUM_THREADS);

    if ( !REF.EndsWith(".fastb") ) build_index_if_needed( REF );

    build_repeat_mask_if_needed( REF, NUM_THREADS );	// depends on global K value gK

    return 0;
}

