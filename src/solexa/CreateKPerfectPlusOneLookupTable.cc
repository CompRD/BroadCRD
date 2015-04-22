#include "MainTools.h"
#include "Basevector.h"
#include "Bitvector.h"

#include <fstream>

ofstream fout;

int main(int argc, char **argv) {
	RunTime();

	BeginCommandArguments;
	CommandArgument_String(KP1TABLE);
	CommandArgument_String(REFERENCE);
	CommandArgument_Int_OrDefault(K, 20);
	CommandArgument_Int_OrDefault(CHUNK, 500000);
	EndCommandArguments;

	vecbitvector bits(KP1TABLE);
	vecbasevector reference(REFERENCE);

	Assert(bits.size() == reference.size());
	Assert(bits[0].size() == reference[0].size());

	basevector kmer;
	int count = 0;
	int filenumber = 0;

	for (size_t i = 0; i < reference.size(); i++) {
		for (unsigned int j = 0; j < reference[i].size(); j++) {
			if (static_cast<int>(bits[i][j]) == 1) {
				if (j+K < reference[i].size()) {
					if (count == 0) {
						fout.open(("ref_" + ToString(filenumber)).c_str());
					}

					kmer.SetToSubOf(reference[i], j, K);
					fout << kmer.ToString() << " " << as_base(reference[i][j+K]) << endl;

					if (count == CHUNK) {
						fout.close();
						count = 0;
						filenumber++;
					} else {
						count++;
					}
				}
			}
		}
	}

	return 0;
}
