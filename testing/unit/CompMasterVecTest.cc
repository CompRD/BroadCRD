//Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

/** testcmv.cc File to test compmastervec.
@class TestCMV
\author: Pablo Alvarez

Things to do for testing:

-fill up a compmastervec with many repeats
-time to see that it is really O(N + KlogK) to fill up.
-test the iterator class and iterator methods (begin, end, both const and non-const).
-test the reference class, in particular that we can use
  cmv[4] = "blah" correctly.

  -test writing and reading one cmv (note that this pretty much requires an operator==).

*/


#include "MainTools.h"
#include "CompMasterVec.h"
#include "CompMasterVecTemplate.h"
#include "CompVecString.h"
#include "testing/TestHelpers.h"
#include "Intvector.h"


/** Testcmv: class for testing compmastervec, with access to cmv internals.

@class Testcmv

\author Pablo Alvarez

This class has been declared a friend of compmastervec so it has access to
private data members. Its methods are used to test different portions of
the compmastervec interface.
*/

struct Testcmv {
  void CreateStrings(int NSTR, int STRSIZE, vecString & strings) {
    strings.Reserve(STRSIZE * NSTR + 10, NSTR);
    String s;
    s.resize(20);
    for (int i = 0; i != NSTR; ++i) {
      StringSetRandom(s, STRSIZE);
      strings.push_back(s);
    }
    cout << "strings.size() is " << strings.size() << endl;
  }    


  double TestTiming(size_t N, ///< total number of items
                    size_t K, ///< number of unique items
                    const vecString & strings) ///< items to insert
  {
    ForceAssertGe(strings.size(), K);
    compvecString cmv;
    cmv.Reserve(0, 0, 0);   
    struct timeval start, end;
    gettimeofday(&start, NULL);
    int j;
    for (size_t i = 0; i != N; ++i) {
      //cycle through one of K unique items to input, avoid using % K.
      j = i%K;
      cmv.push_back(strings[j]);
    }
    gettimeofday(&end, NULL);
    for (size_t i = 0; i != N; ++i) {
      if ( (const String &)(cmv[i]) != strings[i % K]) {
        cout << "comparison failed for i="<<i<<", K="<<K<<endl;
        cout << cmv[i] << ", " << strings[i%K] << endl;
        ForceAssert(0);
      }
    }
    return double(end.tv_sec - start.tv_sec + 
                  (end.tv_usec - start.tv_usec) /1000000.0);
    
  }

  int TestReferenceIterators(const vecString & strings) {
    compvecString cmv;
    const int N = 100;
    const int K = 10;
    const int SIZE = strings[0].size();
    cmv.Reserve(K * SIZE + 10, K, N);
    for (int i = 0; i != N; ++i) {
      cmv.push_back(strings[i % K]);
    }
    compvecString::iterator b = cmv.begin();
    compvecString::iterator e = cmv.end();
    compvecString::iterator i1 = cmv.begin();
    cout << " Test pre-increment " << endl;
    for (int i = 0; i != N; ++i, ++i1) {
      ForceAssertEq((const String &)(*i1), strings[i%K]);
    }
    cout << "Test post-increment " << endl;
    i1 = cmv.begin();
    for (int i = 0; i != N; ++i) {
      ForceAssertEq((const String &)(*i1++), strings[i%K]);
    }
    cout << "Test pre-decrement " << endl;
    i1 = cmv.end();
    for (int i = N-1; i != 0; --i) {
      ForceAssertEq((const String &)(*--i1), strings[i%K]);
    }
    cout << "Test post-decrement " << endl;
    i1 = cmv.end() -1;
    for (int i = N-1; i != 0; --i) {
      ForceAssertEq((const String &)(*i1--), strings[i%K]);
    }
    cout << "Test addition operator " << endl;
    for (int i = 0; i != N; ++i) {
      i1 = cmv.begin() + i;
      ForceAssertEq((const String &)(*i1), strings[i % K]);
    }
    cout << "Test subtraction operator " << endl;
    for (int i = 1; i <= N; ++i) {
      i1 = cmv.end() - i;
      ForceAssertEq((const String &)(*i1), strings[(N-i) % K]);
    }
    cout << "Test reference for assignment " << endl;
    ForceAssert((const String &)cmv[0] != (const String &)cmv[1]);
    //cout << cmv[0] << " " << cmv[1] << endl
    //     << strings[0] << " " << strings[1] << endl;
    cmv[0] = strings[1];
    //cout << cmv[0] << " " << cmv[1] << endl
    //     << strings[0] << " " << strings[1] << endl;
    ForceAssertEq((const String &)cmv[0], strings[1]);
    ForceAssertEq((const String &)cmv[1], strings[1]);

    cout << "All reference and iterator tests passed" << endl;
    return 0;
  }

  int TestAssignmentAppend(const vecString & strings) {
    compvecString cmv;
    const int N = 100;
    const int K = 10;
    const int SIZE = strings[0].size();
    cmv.Reserve(K * SIZE + 10, K, N);
    for (int i = 0; i != N; ++i) {
      cmv.push_back(strings[i % K]);
    }

    cout << "Test assignment and equality operators" << endl;
    compvecString copy1, v3, v4, v5;
    copy1 = cmv;
    ForceAssert(cmv == copy1);

    cout << "Testing Append " << endl;

    int from = 0;
    int to = 5;
    v3.Append(cmv,from, to);
    vec<int> entries;
    for (int i = from; i != to; ++i) {
      entries.push_back(i);
    }
    v4.Append(cmv, entries);
    ForceAssert(v4 == v3);

    v4.Append(cmv, to, cmv.size());
    ForceAssert(v4 == cmv);
      
    entries.clear();
    for (size_t i = 0; i < cmv.size(); i +=2) {
      entries.push_back(i);
    }
    v5.Append(cmv,entries);
    for (size_t i = 0; i < cmv.size(); i +=2) {
      ForceAssert(cmv[i].theX() == v5[i/2].theX());
    }

    /*
    entries.clear();
    int old = v5.size();
    for (int i = 1; i < cmv.size(); i +=2) {
      entries.push_back(i);
    }
    v5.Append(cmv,entries);
    for (int i = 1; i < cmv.size(); i +=2) {
      ForceAssert(cmv[i] == v5[i/2 + old]);
    }

    cout << "All  assignment and append tests passed" << endl;
    */
    return 0;
  }

  int TestFind(const vecString & strings) {
    size_t const NOT_FOUND = static_cast<size_t>(-1L);
    size_t K = 100;
    ForceAssertLe(K, strings.size());
    size_t N = K;
    compvecString cmv;
    const int SIZE = strings[0].size();
    cmv.Reserve(K * SIZE + 10, K, N);
    for (size_t i = 0; i != N; ++i) {
      cmv.push_back(strings[i % K]);
    }
    ForceAssert(cmv.Contains(strings[0]));
    String s;
    StringSetRandom(s,strings[0].size() + 1);
    ForceAssertEq(cmv.Contains(s), false);
    size_t i = cmv.Find(s);
    ForceAssertEq(i, NOT_FOUND);
    i = cmv.UniqueIndex(s);
    ForceAssertEq(i, NOT_FOUND);


    i = cmv.Find(strings[0]);
    ForceAssertEq(i, 0u);
    i = cmv.FindUnique(strings[0]);
    ForceAssertEq(i, 0u);
    i = cmv.UniqueIndex(strings[0]);
    ForceAssertEq(i, 0u);
    ForceAssertEq(i, cmv.UniqueIndex(0));

    cmv.push_back(strings[0]);//create a duplicate!
    cmv.push_back(s);
    size_t pos = cmv.Find(s);
    ForceAssertEq(pos, cmv.size()-1);
    ForceAssertEq(cmv.FindUnique(s), NOT_FOUND); //fail because there are duplicates.
    ForceAssertEq(cmv.UniqueIndex(s), cmv.UniqueSize() - 1);
    ForceAssertEq(cmv.UniqueIndex(s), cmv.UniqueIndex(pos));

    
    //Now add a bunch of repeated strings, and try CreateReverseIndex.
    compvecString cmv2;
    size_t N2 = K * 10;
    cmv2.Reserve(K * SIZE + 10, K, N2);
    for (size_t i = 0; i != N2; ++i) {
      cmv2.push_back(strings[i % K]);
    }
    VecIntVec rindex;
    cmv2.CreateReverseIndex(rindex);
    for (size_t i = 0; i != N2; ++i) {
      size_t u = cmv2.UniqueIndex(cmv2[i]);
      ForceAssertEq(u, cmv2.UniqueIndex(i));
      IntVec const& mates = rindex[u];
      ForceAssertEq(mates.size(), 10u);
      for (IntVec::size_type j = 0; j != mates.size(); ++j) {
        ForceAssert(cmv2[i] == cmv2[mates[j]]);
      }
    }
    
    cout << "Contains, Find, UniqueIndex, CreateReverseIndex tests passed " 
         << endl;
    return 0;
  }

  int TestReadWrite(const vecString & strings) {
    size_t K = 100;
    ForceAssertLe(K, strings.size());
    size_t N = K * 10;
    compvecString cmv;
    unsigned int SIZE = strings[0].size();
    cmv.Reserve(K * SIZE + 10, K, N);
    for (size_t i = 0; i != N; ++i) {
      cmv.push_back(strings[i % K]);
    }
    const String fname="testcmvrw.dat";
    String command = "rm -f " + fname + "*";
    system(command.c_str());

    cmv.WriteAll(fname);
    ForceAssert(IsGoodCompmastervecFile(fname));
    cout << "File written successfully  " << endl;

    compvecString cmv2;
    cmv2.ReadAll(fname);
    ForceAssert(cmv == cmv2);
    cout << "File read successfully  " << endl;

    compvecString cmv3(fname);
    ForceAssert(cmv == cmv3);   
    cout << "Constructor from file successful " << endl;

    String f2 = "partial." + fname;
    int from = 0;
    int to = cmv.size() / 2;
    Remove(f2);
    cmv.Write(f2, from, to);
    compvecString cmv4(f2);
    compvecString cmv5;
    cmv5.Append(cmv, from, to);
    ForceAssert(cmv4 == cmv5);
    cout << "Partial write to file successful " << endl;

    compvecString cmv6;
    cmv6.Append(cmv,0,3);
    cmv6.Append(cmv,0+K, 5+K);
    cout << cmv6 << endl;
    cout << "Output to stdout successful " << endl;

    return 0;
  }        

  int TestPermute(const vecString & strings) {
    size_t K = 10;
    unsigned int SIZE = strings[0].size();
    ForceAssertLe(K, strings.size());
    compvecString from;
    from.Reserve(K * SIZE + 10, K, K);
    for (size_t i = K; i > 0; ) {
      from.push_back(strings[--i]);
    }
    vec<String> to(K/2);
    for (int i = 0; i != (longlong) to.size(); ++i) {
      to[i] = from[i*2];
    }
    to.push_back("notfound");
    //cout << to << "\n\ncompmastervec:\n" << from << endl;
    vec<int> perm;
    vec<int> notFound;
    from.FindPermutationTo(to, perm, notFound);
    cout << "Permutation: " << perm << endl;
    cout << "notFound: " << notFound << endl;
    for (unsigned int i = 0; i != K; ++i) {
      if (perm[i] < (longlong) to.size() && perm[i] != -1) {
        Assert( from[i].theX() == to[perm[i]] );
      }
    }
    PRINT(from);
    PermuteSwappableVec(from, perm);
    cout << "Permuted from: " << endl;
    PRINT(from);
    cout << "All permutation tests passed" << endl;
    return 0;
  }
     

}; //end of class testcmv

int main(int argc, char ** argv) {
  Testcmv t;
  int NSTR = argc > 1 ? atoi(argv[1]): 65000;
  int STRSIZE = argc > 2 ? atoi(argv[2]) : 20;
  vecString strings;
  t.CreateStrings(NSTR, STRSIZE, strings);
  t.TestReferenceIterators(strings);
  t.TestAssignmentAppend(strings);
  t.TestReadWrite(strings);
  t.TestFind(strings);
  t.TestPermute(strings);
  return 0;
}
