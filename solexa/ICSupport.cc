#include "solexa/ICSupport.h"

ICInfo isIC(basevector &read, vecbasevector &fwcontrols, int IC_SEARCH_LENGTH, int MAX_ERR) {
	ICInfo *fwicinfo = new ICInfo[fwcontrols.size()];
	ICInfo *rcicinfo = new ICInfo[fwcontrols.size()];
	ICInfo besticinfo;
	basevector rccontrol;

	for (size_t controlnum = 0; controlnum < fwcontrols.size(); controlnum++) {
		fwicinfo[controlnum].control = controlnum;
		rcicinfo[controlnum].control = controlnum;
		rcicinfo[controlnum].isRC = 1;

		rccontrol = fwcontrols[controlnum];
		rccontrol.ReverseComplement();

		for (int i = 0; i < ((read.isize() < IC_SEARCH_LENGTH) ? read.isize() : IC_SEARCH_LENGTH); i++) {
			if (read[i] != fwcontrols[controlnum][i]) { fwicinfo[controlnum].editDistance++; }
			if (read[i] != rccontrol[i]) { rcicinfo[controlnum].editDistance++; }
		}

		if (fwicinfo[controlnum].editDistance <= MAX_ERR || rcicinfo[controlnum].editDistance <= MAX_ERR) {
			besticinfo = (fwicinfo[controlnum].editDistance < rcicinfo[controlnum].editDistance) ? fwicinfo[controlnum] : rcicinfo[controlnum];
			besticinfo.isIC = 1;
		}
	}

	delete [] fwicinfo;
	delete [] rcicinfo;

	return besticinfo;
}

vec<ICInfo> findICs(vecbasevector &reads, vecbasevector &fwcontrols, int IC_SEARCH_LENGTH, int MAX_ERR) {
	vec<ICInfo> ICs(reads.size());
	int index = 0;
	ICInfo icinfo;

	for (size_t i = 0; i < reads.size(); i++) {
		icinfo = isIC(reads[i], fwcontrols, IC_SEARCH_LENGTH, MAX_ERR);
		if (icinfo.isIC) {
			ICs[index] = icinfo;
			ICs[index].queryId = i;

			index++;
		}
	}

	ICs.resize(index);

	return ICs;
}

vec<short int> getReadMaskFromICs(String icmatchpath, int numReads) {
	vec<short int> readmask(numReads, 0);
	ifstream icmatchfile(icmatchpath.c_str());
	int index = 0;
	short int buffer;

	while (icmatchfile >> buffer) {
		readmask[index] = buffer;
		index++;
	}

	icmatchfile.close();

	return readmask;
}

bitvector getReadBitMaskFromICs(String icmatchpath, int numReads) {
	bitvector readmask(numReads, 0);
	ifstream icmatchfile(icmatchpath.c_str());
	int index = 0, buffer;

	while (icmatchfile >> buffer) {
		readmask.Set(index, ((buffer == -1) ? 0 : 1));
		index++;
	}

	icmatchfile.close();

	return readmask;
}

vec< pair<float, float> > getKPlus1ErrorRate(vecbasevector &reads, vec<ICInfo> &ICs, vecbasevector &fwcontrols, int K) {
	basevector tempcontrol, tempread;
	int smallerlength = (fwcontrols[0].size() < reads[0].size() ? fwcontrols[0].size() : reads[0].size());
	vec< pair<float, float> > errorRate(smallerlength);

	for (unsigned int i = 0; i < errorRate.size(); i++) {
		errorRate[i].first = 0.0;
		errorRate[i].second = 0.0;
	}

	for (unsigned int icindex = 0; icindex < ICs.size(); icindex++) {
		for (int i = 0; i < smallerlength - K; i++) {
			tempcontrol = fwcontrols[ICs[icindex].control];
			if (ICs[icindex].isRC) { tempcontrol.ReverseComplement(); }
			tempcontrol.SetToSubOf(tempcontrol, i, K);

			tempread.SetToSubOf(reads[ICs[icindex].queryId], i, K);

			if (tempcontrol == tempread) {
				base_t baseInQuestion = reads[ICs[icindex].queryId][i + K];
				base_t baseReference = (ICs[icindex].isRC) ? GetComplementaryBase(fwcontrols[ICs[icindex].control][fwcontrols[ICs[icindex].control].size() - i - K - 1]) : fwcontrols[ICs[icindex].control][i + K];

				errorRate[i + K].first += static_cast<float>(baseInQuestion != baseReference);
				errorRate[i + K].second += 1.0;
			}

		}
	}

	return errorRate;
}

vec<float> getTraditionalErrorRate(vecbasevector &reads, vec<ICInfo> &ICs, vecbasevector &fwcontrols) {
	basevector tempcontrol;
	int smallerlength = (fwcontrols[0].size() < reads[0].size() ? fwcontrols[0].size() : reads[0].size());
	vec<float> errorRate(smallerlength, 0.0);

	for (unsigned int icindex = 0; icindex < ICs.size(); icindex++) {
		tempcontrol = fwcontrols[ICs[icindex].control];
		if (ICs[icindex].isRC) { tempcontrol.ReverseComplement(); }

		for (int cycle = 0; cycle < smallerlength; cycle++) {
			errorRate[cycle] += static_cast<float>(reads[ICs[icindex].queryId][cycle] != tempcontrol[cycle]);
		}
	}

	for (int cycle = 0; cycle < smallerlength; cycle++) {
		errorRate[cycle] /= ICs.size();
	}

	return errorRate;
}

