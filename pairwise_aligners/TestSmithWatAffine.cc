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
#include "PrintAlignment.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "FastaFileset.h"
#include <iostream>
#include <fstream>

int main( int argc, char** argv )
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_String_OrDefault(FASTA,"");
    CommandArgument_String_OrDefault(FASTB,"");
    EndCommandArguments;

    vecbvec reads;

    if ( FASTA != "" ) {
	FastFetchReads( reads, nullptr, FASTA );
    } else if ( FASTB != "" ) {
	reads.ReadAll(FASTB);
    } else
	FatalErr("must set either FASTA or FASTB");

    for ( size_t i = 0; i < reads.size()-1; ++i ) {
	for  ( size_t j = i+1; j < reads.size(); ++j ) {
	    cout << "Align " << i << " <- " << j << endl;
	    alignment al;
	    int cost = SmithWatAffineParallel2(reads[i], reads[j], al, true, false);
	    align a = align(al);

	    PrintVisualAlignment(True, cout, reads[i], reads[j], a );
	    cout << "ALIGNMENT: " << endl;
	    cout << al << endl;
	    avector<int> gaps, lens;
	    int pos1, pos2;
	    int errors;
	    al.Unpack( pos1, pos2, errors, gaps, lens );
	    cout << "GAPS LENS" << endl;
	    for ( size_t i = 0; i < gaps.length; ++i )
		cout << gaps(i) << " " << lens(i) << endl;
	    cout << "cost=" << cost;
	    cout << "errors=" << errors << endl;
            cout << endl;
	}
    }


    return 0;
}
