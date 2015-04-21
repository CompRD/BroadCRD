#ifndef _IC_SUPPORT_H
#define _IC_SUPPORT_H

#include "Basevector.h"
#include "Bitvector.h"

class ICInfo {
	public:
		ICInfo() : isIC(0), isRC(0), control(-1), editDistance(0), queryId(0) {}
		bool isIC;
		bool isRC;
		int control; // TODO: potential truncation of index
		int editDistance;
		unsigned long int queryId;
};

ICInfo isIC(basevector &read, vecbasevector &fwcontrols, int IC_SEARCH_LENGTH = 36, int MAX_ERR = 12);
vec<ICInfo> findICs(vecbasevector &reads, vecbasevector &fwcontrols, int IC_SEARCH_LENGTH = 36, int MAX_ERR = 12);
bitvector getReadBitMaskFromICs(String icmatchpath, int numReads);
vec<short int> getReadMaskFromICs(String icmatchpath, int numReads);
vec< pair<float, float> > getKPlus1ErrorRate(vecbasevector &reads, vec<ICInfo> &ICs, vecbasevector &fwcontrols, int K = 20);
vec<float> getTraditionalErrorRate(vecbasevector &reads, vec<ICInfo> &ICs, vecbasevector &fwcontrols);

#endif
