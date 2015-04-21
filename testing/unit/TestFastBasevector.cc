#include "MainTools.h"
#include "FastBasevector.h"
#include <climits>
#include <iostream>

using namespace std;

int main(int argc, char **argv) {
	RunTime();

	//BeginCommandArguments;
	//CommandArgument_String(HEAD);
	//EndCommandArguments;

	cout << "Testing..." << endl;

	fastbasevector a;

	Assert(FAST_BASE_A^1 == FAST_BASE_T);
	Assert(FAST_BASE_T^1 == FAST_BASE_A);
	Assert(FAST_BASE_G^1 == FAST_BASE_C);
	Assert(FAST_BASE_C^1 == FAST_BASE_G);

	//printf("%.4hx %.4hx %.4hx\n", FAST_BASE_A, FAST_BASE_A^1, FAST_BASE_T);
	//printf("%.4hx %.4hx %.4hx\n", FAST_BASE_C, FAST_BASE_C^1, FAST_BASE_G);
	//printf("%.4hx %.4hx %.4hx\n", FAST_BASE_T, FAST_BASE_T^1, FAST_BASE_A);
	//printf("%.4hx %.4hx %.4hx\n", FAST_BASE_G, FAST_BASE_G^1, FAST_BASE_C);

	Assert(fastbasevector::as_base(FAST_BASE_A^1) == 'T');
	Assert(fastbasevector::as_base(FAST_BASE_C^1) == 'G');
	Assert(fastbasevector::as_base(FAST_BASE_T^1) == 'A');
	Assert(fastbasevector::as_base(FAST_BASE_G^1) == 'C');

	Assert(a.size() == 0);
	//Assert(a[0] == FAST_BASE_UNKNOWN);
	//Assert(fastbasevector::as_base(a[0]) == '?');

	String kmer = String("AGTC");
	fastbasevector b(kmer);

	Assert(b.size() == 4);
	Assert(b[0] == FAST_BASE_A);
	Assert(b[1] == FAST_BASE_G);
	Assert(b[2] == FAST_BASE_T);
	Assert(b[3] == FAST_BASE_C);
	//Assert(b[4] == FAST_BASE_UNKNOWN);
	Assert(fastbasevector::as_base(b[0]) == 'A');
	Assert(fastbasevector::as_base(b[1]) == 'G');
	Assert(fastbasevector::as_base(b[2]) == 'T');
	Assert(fastbasevector::as_base(b[3]) == 'C');
	//Assert(fastbasevector::as_base(b[4]) == '?');

	b.ReverseComplement();

	Assert(b.size() == 4);
	Assert(b[0] == FAST_BASE_G);
	Assert(b[1] == FAST_BASE_A);
	Assert(b[2] == FAST_BASE_C);
	Assert(b[3] == FAST_BASE_T);
	//Assert(b[4] == FAST_BASE_UNKNOWN);
	Assert(fastbasevector::as_base(b[0]) == 'G');
	Assert(fastbasevector::as_base(b[1]) == 'A');
	Assert(fastbasevector::as_base(b[2]) == 'C');
	Assert(fastbasevector::as_base(b[3]) == 'T');
	//Assert(fastbasevector::as_base(b[4]) == '?');

	fastbasevector c("AGCTC");
	c.ReverseComplement();
	fastbasevector d("GAGCT");

	for (unsigned int i = 0; i < c.size(); i++) {
		Assert(c[i] == d[i]);
	}

	fastbasevector e("AAGCTA");
	fastbasevector f(6);
	f = e;

	Assert(e[0] == FAST_BASE_A);
	Assert(e[1] == FAST_BASE_A);
	Assert(e[2] == FAST_BASE_G);
	Assert(e[3] == FAST_BASE_C);
	Assert(e[4] == FAST_BASE_T);
	Assert(e[5] == FAST_BASE_A);

	fastbasevector g = e;

	Assert(g[0] == FAST_BASE_A);
	Assert(g[1] == FAST_BASE_A);
	Assert(g[2] == FAST_BASE_G);
	Assert(g[3] == FAST_BASE_C);
	Assert(g[4] == FAST_BASE_T);
	Assert(g[5] == FAST_BASE_A);

	Assert(g == g);
	Assert(g == e);
	Assert(g != b);

	vecbasevector reads("/seq/solexaproc/pipelineOutput/080215_2036V/2036V.4.fastb");
	for (size_t i = 0; i < reads.size(); i++) {
		basevector bv;
		fastbasevector fbv;

		bv.SetToSubOf(reads[i], 0, 20);
		fbv.SetToSubOf(reads[i], 0, 20);

		//bv.ReverseComplement();
		//fbv.ReverseComplement();

		for (int j = 0; j < 20; j++) {
			Assert(as_base(bv[j]) == fastbasevector::as_base(fbv[j]));
		}
	}

	vecbasevector ref1("/seq/solexaproc/pipelineOutput/080215_2036V/2036V.4.fastb");
	vecfastbasevector ref2("/seq/solexaproc/pipelineOutput/080215_2036V/2036V.4.fastb");
	for (size_t i = 0; i < ref1.size(); i++) {
		for (unsigned int j = 0; j < ref1[i].size(); j++) {
			AssertEq(as_base(ref1[i][j]), fastbasevector::as_base(ref2[i][j]));
		}
	}

	cout << "Testing complete!" << endl;

	return 0;
}
