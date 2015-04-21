#include "MainTools.h"
#include "Basevector.h"
#include "solexa/ICSupport.h"

#include <fstream>
#include <map>

using namespace std;

ofstream fout, ferr;

bool fequal(float a, float b) {
	return (fabs(a - b) < 0.0001);
}

int main(int argc, char **argv) {
	RunTime();

	BeginCommandArguments;
	CommandArgument_String(HEAD);
	CommandArgument_Int_OrDefault(K, 20);
	CommandArgument_Int_OrDefault(IC_SEARCH_LENGTH, 36);
	CommandArgument_Int_OrDefault(MAX_ERR, 12);
	CommandArgument_String_OrDefault(OUT, "/dev/stdout");
	CommandArgument_String_OrDefault(ERR, "/dev/stderr");
	CommandArgument_String_OrDefault(CONTROL_REFERENCE_PATH, "/seq/references/Synthetic_internal_controls_set1/v0/Synthetic_internal_controls_set1.fasta.lookuptable.fastb");
	CommandArgument_Bool_OrDefault(COMPUTE_TRADITIONAL_ERROR, 0);
	EndCommandArguments;

	fout.open(OUT.c_str());
	ferr.open(ERR.c_str());

	// Initialize
	vecbasevector reads(HEAD + ".fastb");
	vecbasevector fwcontrols(CONTROL_REFERENCE_PATH);

	// Find the ICs
	vec<ICInfo> ICs = findICs(reads, fwcontrols, IC_SEARCH_LENGTH, MAX_ERR);
	fout << "# ic_count=" << ICs.size() << endl;

	// Get Error Rate
	if (COMPUTE_TRADITIONAL_ERROR) {
		vec<float> errorRate = getTraditionalErrorRate(reads, ICs, fwcontrols);

		for (unsigned int i = K; i < errorRate.size(); i++) {
			fout << i << " " << (fequal(errorRate[i], 0.0) ? 0.0 : errorRate[i]) << endl;
		}
	} else {
		vec< pair<float,float> > errorRate = getKPlus1ErrorRate(reads, ICs, fwcontrols, K);

		for (unsigned int i = K; i < errorRate.size(); i++) {
			fout << i << " " << (fequal(errorRate[i].second, 0.0) ? 0.0 : errorRate[i].first/errorRate[i].second) << " " << errorRate[i].second << endl;
		}
	}

	fout.close();
	ferr.close();

	return 0;
}
