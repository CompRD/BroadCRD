//Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology


#include "math/HoInterval.h"
#include "TestAssert.h"

struct HoIntervalTester {
  static int TestMemberSorted () {
    vec<ho_interval> v;
    for (int i=0; i != 10; ++i) {
      v.push_back(ho_interval(-i*10, -i*10+5));
    }
    sort(v.begin(), v.end());
    for (int i=0; i !=10; ++i) {
      TestAssert(MemberSorted(v,-i*10,5));
      //cout << -i*10 << " " << MemberSorted(v,-i*10,10) << endl;
      TestAssert(!MemberSorted(v,-i*10+5,5));
      //cout << -i*10+5 << " " << MemberSorted(v,-i*10+5,10) << endl;
    }
    return 0;
  }

  int TestContains() {
    HoIntervalWithId i1(0,10,1);
    HoIntervalWithId i2(0,10,2);
    HoIntervalWithId i3(0,10,1);
    HoIntervalWithId i4(1,11,1);
    TestAssert(i1.Contains(i3));
    TestAssert(!(i1.Contains(i2)));
    TestAssert(!(i1.Contains(i4)));
    return 0;
  }
};


int main(int argc, char ** argv) {
  HoIntervalTester t;
  t.TestMemberSorted();
}

