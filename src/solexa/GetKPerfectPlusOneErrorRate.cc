#include "MainTools.h"
#include "String.h"
#include "FastBasevector.h"
#include "MapFastBasevector.h"
#include "Bitvector.h"
#include "solexa/KPerfectPlusOneUtils.h"

#include <fstream>

#include <ctime>
#include <cstdlib>
#include <limits>
#include <unistd.h>

#ifdef __GNUC__
#include <ext/hash_map>
using __gnu_cxx::hash_map;
#endif

using namespace std;

class ErrorInfo {
	public:
		ErrorInfo() : m_errorCount(0), m_errorTotal(0) {};

		unsigned int m_errorCount;
		unsigned int m_errorTotal;
};

typedef struct {
	fast_base_t nextbase;
	fast_base_t prevbase;
} boundary_bases;

ofstream fout, ferr;

typedef __gnu_cxx::hash_map<fastbasevector, boundary_bases, fastbasevector_hash, fastbasevector_equal> Map;

Map getKmerMap(String UNIQUE_KMER_LOOKUP_FILE, String REFERENCE, int K,
                size_t CONTIG_START, size_t CONTIG_END,
                unsigned int BASE_START, unsigned int BASE_END ) {
	String ktablePath = UNIQUE_KMER_LOOKUP_FILE;
	if (!IsRegularFile(ktablePath)) {
		ferr << "Error: K" << K << " lookup file '" << ktablePath << "' does not exist or is not a regular file." << endl;
		exit(1);
	}

	vecbitvector bits(ktablePath);
	vecfastbasevector reference(REFERENCE);

	Map kmerMap;
	fastbasevector kmer;
	unsigned long int uniqueKmers = 0;
	int matches = 0;

	if ( CONTIG_START > bits.size() ) CONTIG_START = bits.size();
	if ( CONTIG_END > bits.size() ) CONTIG_END = bits.size();
	for (size_t i = CONTIG_START; i < CONTIG_END; ++i) {
                unsigned int start = BASE_START < bits[i].size() ? BASE_START : bits[i].size();
                unsigned int end = BASE_END < bits[i].size() ? BASE_END : bits[i].size();
		for (unsigned int j = start; j < end; ++j) {
			if (bits[i][j] == 1 && j+K < reference[i].size()) {
				kmer.SetToSubOf(reference[i], j, K);

				boundary_bases b;
				b.nextbase = (j+K < bits[i].size()) ? reference[i][j+K] : FAST_BASE_UNKNOWN;
				b.prevbase = (j > 0) ? reference[i][j-1] : FAST_BASE_UNKNOWN;
				kmerMap[kmer] = b;

				++uniqueKmers;
			}
		}
	}
	fout << "# unique_kmers=" << uniqueKmers << endl;

	return kmerMap;
}

int main(int argc, char **argv) {
	RunTime();

	BeginCommandArguments;
	CommandArgument_String(READS);
    CommandArgument_String(REFERENCE);
	CommandArgument_String(UNIQUE_KMER_LOOKUP_FILE);
	CommandArgument_Int_OrDefault(K, 20);
	CommandArgument_String_OrDefault(OUT, "/dev/stdout");
	CommandArgument_String_OrDefault(ERR, "/dev/stderr");
	//CommandArgument_String_OrDefault(HASHED_KMER_LOOKUP_PATH, "/home/radon01/kiran/tmp/kp1lookup");
	CommandArgument_UnsignedInt_OrDefault(CONTIG_START, 0);
	CommandArgument_UnsignedInt_OrDefault(CONTIG_END, INT_MAX);
	CommandArgument_UnsignedInt_OrDefault(BASE_START, 0);
	CommandArgument_UnsignedInt_OrDefault(BASE_END, INT_MAX);
	CommandArgument_Int_OrDefault(REF_PIECE_START, 0);
	CommandArgument_Int_OrDefault(REF_PIECE_END, 1);
	CommandArgument_Bool_OrDefault(CHUNK, 0);
	EndCommandArguments;

	fout.open(OUT.c_str());
	ferr.open(ERR.c_str());

	// Initialize
	vecfastbasevector reads(READS);
	fastbasevector kmer(K);
	fast_base_t nextbase;

	vec<ErrorInfo> errorRate(reads[0].size());
	unsigned int matches = 0;
	unsigned int numReads = reads.size();
	unsigned int readLength = reads[0].size();
	unsigned int numKmers = reads[0].size() - K;
	unsigned int jplusk;

	// Processing kmers
	for (int chunk = REF_PIECE_START; chunk < REF_PIECE_END; chunk++) {
		Map kmerMap = getKmerMap(UNIQUE_KMER_LOOKUP_FILE, REFERENCE, K, CONTIG_START, CONTIG_END, BASE_START, BASE_END);

		for (unsigned int i = 0; i < numReads; ++i) {
			for (unsigned int j = 0; j < numKmers; ++j) {
				kmer.SetToSubOf(reads[i], j, K);
				jplusk = (j+K);

				if (kmerMap.count(kmer) > 0) {
					++matches;

					errorRate[jplusk].m_errorCount += (reads[i][jplusk] != kmerMap[kmer].nextbase);
					++errorRate[jplusk].m_errorTotal;
				}
			}

			reads[i].ReverseComplement();

			for (unsigned int j = 1; j <= numKmers; ++j) {
				kmer.SetToSubOf(reads[i], j, K);
				
				if (kmerMap.count(kmer) > 0) {
					++matches;

					errorRate[readLength - j].m_errorCount += (reads[i][j-1] != kmerMap[kmer].prevbase);
					++errorRate[readLength - j].m_errorTotal;
				}
			}
		}
	}

	fout << "# matches=" << matches << endl;

	for (unsigned int i = 0; i < errorRate.size(); ++i) {
		fout << i << " " << errorRate[i].m_errorCount << " " << errorRate[i].m_errorTotal << endl;
	}

	fout.close();
	ferr.close();

	return 0;
}
