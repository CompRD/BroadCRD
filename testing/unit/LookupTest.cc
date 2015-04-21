/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "lookup/Hit.h"
#include "TestAssert.h"
#include "random/Random.h"

class LookupTester {
  int verbose;

public:

  LookupTester(int verbose): verbose(verbose) { }
  
  void TestBasesToQueries() {
    const unsigned int K = 12;
    basevector b(2*K);
    for ( unsigned int i=0; i<2*K; ++i )
      b.Set(i,0);
    b.Set(K,1);
    vec<Query> r;
    BasesToQueries(b, r, K);
    TestAssertEq(r.size(), b.size()-(K-1));
    TestAssertEq(r[0].Kmer(), 0u);
    TestAssertEq(r[0].QueryPos(), 0);
    unsigned int exp=1;
    for (unsigned int i=1; i<r.size(); ++i, exp*=4) {
      TestAssertEq(r[i].QueryPos(), int(i));
      TestAssertEq(r[i].Kmer(), exp);
    }
    cout << "BasesToQueries test passed." << endl;
  }

  void TimeTestBasesToQueries() {
    const unsigned int N = 1000000;
    const unsigned int K = 12;
    basevector b(N);
    for ( unsigned int i = 0; i < N; i++ )
      b.Set( i, randomx( ) % 4 );
    vec<Query> r;
    BasesToQueries(b, r, K);
    r.clear();
    double startTime = WallClockTime();
    for ( unsigned int i=0; i< N - (K-1); ++i )
      r.push_back(Query(i, Index(b, i, K)));
    cout << "Simple generation of query sequence took "
	 << TimeSince(startTime) << endl;
    r.clear();
    startTime = WallClockTime();
    BasesToQueries(b, r, K);
    cout << "BasesToQueries took " << TimeSince(startTime) << endl;
  }

  void TestVecBasesToQueries() {
    const unsigned int K = 12;
    basevector b(2*K);
    for ( unsigned int i=0; i<2*K; ++i )
      b.Set(i,0);
    b.Set(K,1);
    
    vecbasevector v;
    const unsigned int N = 10;
    for (unsigned int i=0; i<N; ++i)
      v.push_back(b);
    
    vec<Query> r;
    BasesToQueries(v, r, K);
    unsigned int Q = b.size()-(K-1);
    TestAssertEq(r.size(), N*Q);
    TestAssertEq(r[0].Kmer(), 0u);
    TestAssertEq(r[0].QueryPos(), 0);
    for (unsigned int s=0; s<N; ++s) {
      unsigned int exp=1;
      for (unsigned int i=1; i<Q; ++i, exp*=4) {
	TestAssertEq(r[s*Q+i].QueryPos(), int(s*b.size()+i));
	TestAssertEq(r[s*Q+i].Kmer(), exp);
      }
    }
    cout << "VecBasesToQueries test passed." << endl;
  }

};

int main(int argc, char ** argv) {
  cout << "\nLookupTest:\n";
  LookupTester t(argc-1);
  t.TestBasesToQueries();
  t.TestVecBasesToQueries();
  t.TimeTestBasesToQueries();
  return 0;
}
