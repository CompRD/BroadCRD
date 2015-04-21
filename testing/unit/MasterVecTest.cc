/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
// MakeDepend: archived

#include "MainTools.h"
#include "VecString.h"
#include "VecUtilities.h"
#include "TaskTimer.h"
#include "testing/TestHelpers.h"
#include "TestAssert.h"
#include "Intvector.h"
#include "random/Random.h"

#include <fcntl.h>
#include <iostream>
#include <cmath>

using namespace std;


/** Testmv. Class for some simple tests of mastervec function.
 * @class Testmv.
 * \author Pablo Alvarez
 * See also VecStringTest.cc
 */

class Testmv {
 public:
  int TestPushBackReserve() {
    const int K = 4;
    vecString vecs[K];
    const int N = 10000000;
    vecString strings;
    strings.Reserve(N * 21, N);
    String s;
    for (int i = 0; i != N; ++i) {
      StringSetRandom(s, 20);
      strings.push_back(s);
    }
    int max = N;
    TaskTimer t;
    vec<float> times(K);
    const double FACTOR = 4.0;
    for (int i = 0; i != K; ++i) {
      t.Start();
      for (int j = 0; j != max; ++j) {
        vecs[i].push_back_reserve(strings[j]);
      }
      t.Stop();
      times[i] = t.Elapsed();
      t.Reset();
      max = int(max/ FACTOR);
    }
    cout << "data size:";
    for (int i = 0; i != K; ++i) cout << vecs[i].size() << "  ";
    cout << "\ntime:   ";
    copy(times.begin(), times.end(), ostream_iterator<float>(cout, "  "));
    //Verify that insertion has been between O(N) and O(NlogN).
    for (int i = 0; i != K -1; ++i) {
      double top = 1.3 * (times[i] / FACTOR);
      double bottom = 0.7 * (times [i] / (FACTOR * log(FACTOR)));
      if (bottom < 0.005) break; //too small, not enough accuracy.
      Assert(top > times[i+1]);
      Assert(bottom < times[i+1]);
    }
    cout << "\npush_back_reserve tests passed." << endl;
    return 0;
  }

  int TestSelfPushBackReserve () {
    VecIntVec v;
    for (int i=0; i != 10000; ++i) {
      int s = randomx() % 10;
      IntVec serf;
      for (int j=0; j != s; ++j) serf.push_back(randomx());
      v.push_back_reserve(serf);
    }
    VecIntVec first(v), second(v);
    for (int i=0; i != 1000000; ++i) {
      DotMod(cout, i, 100000);
      int pos = randomx() % v.size();
      first.push_back_reserve(v[pos]);
      v.push_back_reserve(v[pos]);
      second.push_back_reserve(v[pos]);
    }
    for (VecIntVec::size_type i=0; i != v.size(); ++i) {
      if (v[i] != first[i] || v[i] != second[i]) {
	//PRINT4(i,v[i], first[i], second[i]);
	TestAssert(v[i] == first[i]);
	TestAssert( v[i] == second[i]);
      }
    }
    cout << "push_back_reserve own element tests passed" << endl;
    return 0;
  }

  int TestReadWrite() {
    vecString a, b;
    String fname = "test.dat";
    a.push_back("String 1");
    a.push_back("String 2");
    a.push_back("String 3");
    a.WriteAll(fname);
    b.ReadAll(fname);
    cout << " raw count " << MastervecFileRawCount(fname, sizeof(char)) << endl;
    ForceAssertEq(MastervecFileObjectCount(fname), 3);
    ForceAssert(a == b);
    cout << a << endl;
    cout << "ReadWrite tests passed " << endl;
    return 0;
  }

  int TestRawsize() {
    vec<String> strings(10);
    vecString v1, v2, v3, v4, v5;
    ForceAssert(v1.empty());
    for (int i = 0; i != (longlong) strings.size(); ++i) {
      StringSetRandom(strings[i],20);
      v1.push_back(strings[i]);
    }
    ForceAssert(!v1.empty());
    v2 = v1;
    vec<int> entries(3);
    entries[0]=0;
    entries[1]=1;
    entries[2]=2;
    vec<int> e2(entries);
    e2[0]=3;

    //Compare rawsize(from,to) to rawsize(entries).
    longlong r1 = v1.rawsize(0,3);
    longlong r2 = v1.rawsize(entries);
    ForceAssertEq(r1,r2);

    //Make sure rawsize(from, to) still works after a swap.
    v2.SwapElements(0,3);
    r1 = v1.rawsize(e2);
    r2 = v2.rawsize(0,3);
    ForceAssertEq(r1,r2);
    cout << "rawsize tests passed " << endl;
    return 0;
  }


  int TestAssignAppend() {
    vec<String> strings(10);
    vecString v1, v2, v3, v4, v5;
    for (int i = 0; i != (longlong) strings.size(); ++i) {
      StringSetRandom(strings[i],20);
      v1.push_back(strings[i]);
    }
    v2 = v1;
    ForceAssert(v2 == v1);
    int from = 0;
    int to = 5;
    v3.Append(v1,from, to);
    vec<int> entries;
    for (int i = from; i != to; ++i) {
      entries.push_back(i);
    }
    v4.Append(v1, entries);
    ForceAssert(v4 == v3);

    //Append from to onto existing mvec
    v4.Append(v1, to, v1.size());
    ForceAssert(v4 == v1);

    //Append entries onto blank mvec
    entries.clear();
    for (int i = 0; i < v1.size(); i +=2) {
      entries.push_back(i);
    }
    v5.Append(v1,entries);
    for (int i = 0; i < v1.size(); i +=2) {
      ForceAssert(v1[i] == v5[i/2]);
    }

    //Append entries onto existing mvec
    entries.clear();
    for (int i = 1; i < v1.size(); i +=2) {
      entries.push_back(i);
    }
    int old = v5.size();
    v5.Append(v1,entries);
    for (int i = 1; i < v1.size(); i +=2) {
      ForceAssert(v1[i] == v5[i/2 + old]);
    }


    cout << "AssignAppend tests passed " << endl;
    return 0;
  }

  int TestPermute() {
    int a[] = {0,1,2,3,4,5,6,7,8,9};
    int p[] = {0,3,1,2,5,4,7,8,9,6};
    vec<int> av(10);copy(a,a+10,av.begin());
    vec<int> pv(10);copy(p,p+10,pv.begin());
    vecString orig, permuted;
    orig.Reserve(av.size() * 10, av.size());
    for (int i = 0; i != (longlong) av.size(); ++i) {
      orig.push_back(ToString(av[i]));
      permuted.push_back(ToString(av[i]));
    }
    PermuteSwappableVec(permuted, pv);
    PermuteVec(av, pv);
    for (int i = 0; i != (longlong) av.size(); ++i) {
      ForceAssertEq(av[i], permuted[i].Int());
      ForceAssertEq(permuted[pv[i]], orig[i]);
    }
    cout << "Permutation tests passed " << endl;
    return 0;
  }

};


int main(int argc, char ** argv) {
  Testmv t;
  t.TestPushBackReserve();
  t.TestSelfPushBackReserve();
  t.TestReadWrite();
  t.TestRawsize();
  t.TestAssignAppend();
  return t.TestPermute();
}

