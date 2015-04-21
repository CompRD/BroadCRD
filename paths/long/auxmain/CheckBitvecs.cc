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



int main(int argc, char* argv[] )
{
    RunTime( );

    BeginCommandArguments;
    CommandArgument_String(BV1);
    CommandArgument_String(BV2);
    EndCommandArguments;

    vecbitvector bv1(BV1);
    vecbitvector bv2(BV2);

    cout << "bv1.size()=" << bv1.size() << endl;
    cout << "bv2.size()=" << bv2.size() << endl;

//    ForceAssertEq(bv1.size(), bv2.size());
    ForceAssertLe(bv1.size(), bv2.size());
    if ( bv1.size() != bv2.size() )
	cout << "warning: bv1.size()=" << bv1.size() << ", bv2.size()=" << bv2.size() << endl;

    for ( size_t g = 0; g < bv1.size(); ++g )
	ForceAssertEq(bv1[g].size(), bv2[g].size() );

    for ( size_t g = 0; g < bv1.size(); ++g ) {
	cout << "contig " << g << endl;
	for ( size_t i = 0; i < bv1[g].size(); ++i ) {
	    if ( bv1[g][i] != bv2[g][i] )
		cout << g << ":" << i << " " << int(bv1[g][i]) << " " << int(bv2[g][i]) << endl;
	}
    }

    return 0;
}

