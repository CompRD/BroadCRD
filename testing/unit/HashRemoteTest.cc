#include "HashLinearRemoteStorage.h"
#include "Vec.h"
#include "Vec.h"
#include "String.h"
#include "system/Assert.h"
#include "VecString.h"

#include <functional>

class TestHashRemote {
 public:
  struct Hash : public binary_function<String const&, int, int> {
    int operator()(String const& v, int M) {
      if (v.size() == 0) return 0;
      return v[0] % M;
    }
  };

  int Test(int = 0 , char ** = 0) {
    vecString mv;
    HashLinearRemoteStorage<String, Hash> h(&mv, 7);
    const int VECS = 100;
    const int MULTIPLIER = 20;
    const int N = 20;
    vec<int> positions(VECS * MULTIPLIER);
    vec<String> copies(VECS);
    for (int i = 0; i != VECS; ++i) {
      String s(N);
      for (int j = 0; j != N; ++j) {
        s[j] = rand();
      }
      //save a separate copy for checking
      copies[i] = s;
      //insert multiple times to make sure that the searching works correctly.
      for (int j = 0; j != MULTIPLIER; ++j) {
        positions[i + VECS * j] = h.Insert(s);
      }
    }

    cout << "positions: " << positions << endl;

    for (int i = 0; i != VECS; ++i) {
      int j = positions[i];
      if (!(copies[i] == h.GetItem(j))) {
        cout << "table access failed"
             << " copies[" << i << "] != h.GetItem("<<j<<")"<<endl;
        cout << copies[i] <<endl;
        cout << h.GetItem(j) << endl;
        ForceAssert(0);
      }
    }
    cout << "Verified that table access gives back the correct vector" << endl;

    cout << "Insert new succeeded" << endl;
    int mustbe4 = h.Find(mv[4]);  ForceAssertEq(mustbe4, 4);
    cout << "Find existing succeeded" << endl;

    int mustbe5 = h.Insert(mv[5]); ForceAssertEq(mustbe5, 5);
    cout << "Insert existing succeeded" << endl;

    String s(3);
    size_t mustbeEMPTY = h.Find(s);
    size_t i = HashLinearRemoteStorage<String, char, Hash>::EMPTY;
    ForceAssertEq(mustbeEMPTY, i);
    cout << "Find nonexistent failed successfully" << endl;


    //copy(mv.begin(), mv.end(), ostream_iterator<IntVec>(cout, "\n\n"));

    mv.WriteAll("testhash.mvec.dat");

    //Now let's see if rehashing messes anything up

  
  
    return 0;
  }
}; //end of class TestHashRemote

int main(int argc, char ** argv) {
  TestHashRemote t;
  return t.Test(argc, argv);
}
