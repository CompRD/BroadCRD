/////////////////////////////////////////////////////////////////////////////
////                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
////       This software and its documentation are copyright (2008) by the   //
////   Broad Institute/Massachusetts Institute of Technology.  All rights    //
////   are reserved.  This software is supplied without any warranty or      //
////   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
////   can be responsible for its use, misuse, or functionality.             //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "lookup/LookAlign.h"
#include "lookup/PairFinderTools.h"
#include "lookup/SerialQltout.h"

inline void ReserveAndIncrement(vec<float> &distrib, unsigned int value) {
    if (distrib.size() <= value) {
        distrib.resize(static_cast<int>(1.3*static_cast<float>(value)));
    }

    distrib[value] += 1.0;
}

unsigned int GetNormalization(vec<float> &distrib1, vec<float> &distrib2, vec<float> &distrib3) {
    unsigned int sum = 0;

    for (unsigned int i = 0; i < distrib1.size(); i++) { sum += distrib1[i]; }
    for (unsigned int i = 0; i < distrib2.size(); i++) { sum += distrib2[i]; }
    for (unsigned int i = 0; i < distrib3.size(); i++) { sum += distrib3[i]; }

    return sum;
}

void PrintDistribution(vec<float> &distrib, unsigned int sum, String filename) {
    float fsum = static_cast<float>(sum);

    ofstream dout(filename.c_str());
    for (unsigned int i = 0; i < distrib.size(); i++) {
        float value = distrib[i]/fsum;
        if (distrib[i] > 0.0) {
            dout << i << " " << value << endl;
        }
    }
    dout.close();
}

int main(int argc, char **argv) {
    RunTime();

    BeginCommandArguments;
    CommandArgument_String_Abbr_Doc(ALIGNS_END_1,A1,"Alignments for the first ends of the pairs");
    CommandArgument_String_Abbr_Doc(ALIGNS_END_2,A2,"Alignments for the second ends of the pairs");
    CommandArgument_String(OUT_HEAD);
    EndCommandArguments;

    SerialQltout aligns_end_1(ALIGNS_END_1);
    SerialQltout aligns_end_2(ALIGNS_END_2);

    vec<float> innerd(10000);
    vec<float> outerd(10000);
    vec<float> otherd(10000);
    unsigned int innermax = 0, outermax = 0, othermax = 0;
    unsigned int innercount = 0, outercount = 0, othercount = 0;
    unsigned int intercount = 0;

    vec<look_align> aligns1;
    vec<look_align> aligns2;

    bool hasMore = aligns_end_1.NextSet(aligns1) && aligns_end_2.NextSet(aligns2) ;
    // Get the fragment size distribution for unique pairs
    cout << "Computing fragment size distributions for unique pairs... " << flush;
    while( hasMore ) {

        if ( aligns1[0].QueryId() < aligns2[0].QueryId() ) {
	  hasMore = aligns_end_1.NextSet(aligns1); // get next read id from aligns_end_1
	  continue;
	}
        if ( aligns1[0].QueryId() > aligns2[0].QueryId() ) {
	  hasMore = aligns_end_2.NextSet(aligns2); // get next read id from aligns_end_2
	  continue;
	}
	// if we got here we know that query ids in aligns1 and aligns2 are the same 

        // if any of the two paired end alignments is not unique, continue to the next query (read pair) 
        if (aligns1.size() != 1 || aligns2.size() != 1) {
	  hasMore = aligns_end_1.NextSet(aligns1) && aligns_end_2.NextSet(aligns2) ;
	  continue;
	}

        const look_align & la1 = aligns1[0];
	const look_align & la2 = aligns2[0];

	if (PairPlacement(la1, la2) != INTRA) { 
	    // alignments are on different contigs, abort and go get next pair of unique end alignments!
	    ++intercount;
	    hasMore = aligns_end_1.NextSet(aligns1) && aligns_end_2.NextSet(aligns2) ;
	    continue;
	}

	unsigned int fragmentsize = FragmentSize(la1, la2);

	switch (PairDirection(la1, la2)) {
	case INNER:
	    ReserveAndIncrement(innerd, fragmentsize);
	    ++innercount;
	    if (fragmentsize > innermax) { fragmentsize = innermax; }
	    break;
	case OUTER:
	    ReserveAndIncrement(outerd, fragmentsize);
	    ++outercount;
	    if (fragmentsize > outermax) { fragmentsize = outermax; }
	    break;
	default:
	    ReserveAndIncrement(otherd, fragmentsize);
	    ++othercount;
	    if (fragmentsize > othermax) { fragmentsize = othermax; }
	    break;
	};
	hasMore = aligns_end_1.NextSet(aligns1) && aligns_end_2.NextSet(aligns2) ;
            
    }

    unsigned int norm = GetNormalization(innerd, outerd, otherd);

    PrintDistribution(innerd, norm, OUT_HEAD + ".inner.distribution");
    PrintDistribution(outerd, norm, OUT_HEAD + ".outer.distribution");
    PrintDistribution(otherd, norm, OUT_HEAD + ".other.distribution");
    
    cout << "done" << endl;

    PrintWithSep(cout,innercount+outercount+othercount); cout << " unique intra-contig pairs\n";
    cout << ' '; PrintWithSep(cout,innercount); cout << " unique inner pairs\n"; 
    cout << ' '; PrintWithSep(cout,outercount); cout << " unique outer pairs\n";
    cout << ' '; PrintWithSep(cout,othercount); cout << " unique left/right pairs\n";
    PrintWithSep(cout, intercount); cout << " unique inter-contig pairs\n";

    return 0;
}
