///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - April 17, 2013
//

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "paths/long/FriendAlignFinderNaif.h"


int main(int argc, char* argv[])
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_String_Doc(FRIENDS, ".friends file to read");
    CommandArgument_String_Doc(OUT_HEAD, "output head for files");
    CommandArgument_UnsignedInt_OrDefault_Doc(MAX_FREQ,32767,"maximum frequency to allow");
    CommandArgument_UnsignedInt_OrDefault_Doc(THREADS,0,"number of threads to use");
    EndCommandArguments;

    THREADS=configNumThreads(THREADS);
    omp_set_num_threads(THREADS);

    VirtualMasterVec<SerfVec<Friend> > friends( FRIENDS );

    const size_t max_t = 32768;

    cout << "number of records is " << friends.size() << endl;


    vector<size_t> counts(max_t+1, 0U);

    for ( size_t i = 0; i < friends.size(); ++i ) {
	if ( 0 == i % 10000 )
	    cout << i << endl;

	size_t c = friends[i].size();
	if ( c > max_t ) {
	    cout << "warning: value of " << c << " being clamped to " << max_t  << endl;
	    c = max_t;
	}
	counts[c] += 1;
    }

    String output = OUT_HEAD + "_counts.txt";
    ofstream out( output );

    for ( size_t j = 0; j < std::min( size_t(MAX_FREQ)+1, counts.size() ) ; j++ )
	out << j << " " << counts[j] << endl;

    out.close();

    String output2 = OUT_HEAD+"_cumsum.txt";
    ofstream out2( output2 );

    String output3 = OUT_HEAD+"_weighted_cumsum.txt";
    ofstream out3( output3 );

    size_t sum = 0;
    size_t ffsum = 0;
    long double sum_max = friends.size();
    for ( size_t j = 0; j < counts.size(); j++ ) {
	sum += counts[j];
	ffsum += j * counts[j];
	out2 << j << " " << sum / sum_max << endl;
	out3 << j << " " << ffsum << endl;
    }

    out2.close();
    out3.close();


    return 0;

}

