//Copyright (c) 2000-2004 Broad Institute

#include "Map.h"
#include "String.h"
#include "testing/TestHelpers.h"

//#include <iterator>
//#include <ext/hash_set>

/** MapTest.cc Test helper functions in Map.h
\file MapTest.cc
*/

struct StringPtrCompare: public binary_function<const String *, const String *, bool> {
 public:
  bool operator()(const String * lhs,  const String * rhs) {
    return *lhs < *rhs;
  }
};

int main(int argc, char ** argv) {
  std::map<String, int> forward;
  std::map<int, String> reverse;
  vec<String> strings;
  String s;
  s.resize(20);
  for (int i = 0; i != 10; ++i) {
    StringSetRandom(s, 20);
    strings.push_back(s);
    forward[s] = i;
    reverse[i] = s;
  }
  vec<String> fs = keys(forward);
  vec<String> rs = values(reverse);

  vec<const String *> cfs;
  keyPointers(forward, cfs);
  vec<const String *> crs;
  valuePointers(reverse, crs);


  sort(fs.begin(), fs.end());
  sort(rs.begin(), rs.end());
  sort(cfs.begin(), cfs.end(), StringPtrCompare());
  sort(crs.begin(), crs.end(), StringPtrCompare());
  sort(strings.begin(), strings.end());

  for (int i = 0; i != 10; ++i) {
    ForceAssertEq(fs[i], strings[i]);
    ForceAssertEq(rs[i], strings[i]);
    ForceAssertEq(*cfs[i], strings[i]);
    ForceAssertEq(*crs[i], strings[i]);
  }
  cout << "All map helper tests succeeded" << endl;
  return 0;
}

/*
$Log: not supported by cvs2svn $
Revision 1.1  2004/11/04 22:16:03  palvarez
Added or moved unit test files into this directory.

Revision 1.1  2004/11/01 21:06:51  palvarez
Tests for the methods in Map.h
  
*/  
    
