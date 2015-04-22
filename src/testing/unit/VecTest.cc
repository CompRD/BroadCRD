#include "Vec.h"
#include "VecUtilities.h"
#include "STLExtensions.h"
#include "feudal/BinaryStream.h"

/** VecTest.cc File to test vec-related methods
\file VecTest.cc
*/

class VecTest {
 public:
  static int TestPermute() {
    int a[]         = {0,1,2,3,4,5,6,7,8,9};
    int p[]         = {0,3,1,2,5,4,7,8,9,6};
    int apermuted[] = {0,2,3,1,5,4,9,6,7,8};//by hand
    vec<int> av(10);copy(a,a+10,av.begin());
    vec<int> pv(10);copy(p,p+10,pv.begin());
    vec<int> avpermuted(av);
    PermuteVec(avpermuted, pv);
    for (int i = 0; i != 10; ++i) {
      ForceAssertEq(avpermuted[i], apermuted[i]);//check by hand
      ForceAssertEq(avpermuted[p[i]], a[i]);//check theoretical equality
    }
    cout << "Permutation tests passed " << endl;
    return 0;
  }
};

int main(int argc, char ** argv) {
  VecTest t;
  int ret = 0;
  ret |= t.TestPermute();
  return ret;
}
