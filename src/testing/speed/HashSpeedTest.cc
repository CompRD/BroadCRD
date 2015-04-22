/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


///Simple speed test for hash_set, HashSimple and a sorted uniqued vec.
/// \file HashSpeedTest.cc

#include "MainTools.h"
#include "TaskTimer.h"
#include "HashSimple.h"
#include "random/Random.h"
#include <ext/hash_set>

using __gnu_cxx::hash_set;

int main( int argc, char *argv[])
{
  RunTime();
  __gnu_cxx::hash<int> h;
  TaskTimer t;
  for (int size=100000; size <= 100000000; size *=10) {
    vec<int> v(size);
    vec<int> v2;
    //allocate memory here so we are treating the custom hash table 
    //and the vector equally.
    v2.reserve(size);
    HashSimple table(size);
    hash_set<int> table2(size);
    int x, y;

    for (int j=0; j != size; ++j) v[j] = randomx();
    t.Reset();
    t.Start();
    for (int j=0; j != size; ++j) { x = h(v[j]); y= ++x; }
    t.Stop();
    cout << size << " hashes " << t << endl;

    t.Reset();
    t.Start();

    for (int j=0; j != size; ++j) table.Insert(v[j]);
    for (int j=0; j != size; ++j) table.Insert(v[j]);
    for (HashSimple::iterator i=table.begin(); i != table.end(); ++i) {
      //PRINT4(i.pos, i.SIZE, table.data[i.pos], table.capacity());
      x = *i;
    }
    t.Stop();
    cout << size << " hashed linear insert & traverse " << t << endl;
    //#if 0
    t.Reset();
    t.Start();
    for (int j=0; j != size; ++j) table2.insert(v[j]);
    for (int j=0; j != size; ++j) table2.insert(v[j]);
    for (hash_set<int>::iterator i=table2.begin(); i != table2.end(); ++i) {
      x = *i;
    }
    t.Stop();
    cout << size << " hash_set single size insert & traverse " << t << endl;

    t.Reset();
    t.Start();
    table2.clear();
    t.Stop();
    cout << size << " hash_set clearing " << t << endl;

    t.Reset();
    t.Start();
    table2.resize(size*2);
    t.Stop();
    cout << size << " hash_set resizing " << t << endl;

    t.Reset();
    t.Start();
    for (int j=0; j != size; ++j) table2.insert(v[j]);
    for (int j=0; j != size; ++j) table2.insert(v[j]);
    for (hash_set<int>::iterator i=table2.begin(); i != table2.end(); ++i) {
      x = *i;
    }
    t.Stop();
    cout << size << " hash_set double size insert & traverse " << t << endl;
    
    //#endif
    t.Reset();
    t.Start();
    for (int i=0; i != size; ++i) v2.push_back(v[i]);
    for (int i=0; i != size; ++i) v2.push_back(v[i]);
    UniqueSort(v2);
    for (int i=0; i != v2.isize(); ++i) x = v2[i];
    t.Stop();
    cout << size << " vec push_back, UniqueSort & traverse " << t << endl;

    t.Reset();
    t.Start();
    for (int i=0; i != size; ++i) 
      table.Has(v[i]);
    t.Stop();
    cout << size << " hash linear find " << t << endl;

    t.Reset();
    t.Start();
    for (int i=0; i != size; ++i) 
      table2.find(v[i]);
    t.Stop();
    cout << size << " hash_set double size find " << t << endl;

    t.Reset();
    t.Start();
    for (int i=0; i != size; ++i) 
      binary_search(v2.begin(), v2.end(), v[i]);
    t.Stop();
    cout << size << " vec find " << t << endl;
  }
  return 0;
}

