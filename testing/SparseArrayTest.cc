#include "MainTools.h"
#include "TaskTimer.h"

#include "math/SparseArray.h"


typedef SparseArray<int,map<unsigned int, int> > Array;

int main() {

  Array a(10000,0);

  cerr << "sparse array of " << a.size() << " elements;" << endl;
  cerr << "takes up " << a.MemoryUsage() << " bytes;\nStarting random access test:"<< endl;

  a.Set(2,2);
  a.Set(4,4);
  a.Set(8,8);
  a.Set(9,9);
  a.Set(256,2);
  a.Set(348,-4);
  a.Set(865,12);
  a.Set(1321,-25);


  cerr << "2=" << a.Get(2) << " 4=" << a.Get(4) 
       << " 8=" << a.Get(8) << " 9=" << a.Get(9) << endl;
  
  cerr << "2=" << a.Get(256) << " -4=" << a.Get(348) 
       << " 12=" << a.Get(865) << " -25=" << a.Get(1321) << endl;
  
  cerr << "Breaking the sortedness: " << endl;

  
  a.Set(18,64);
  cerr << "64=" << a.Get(18) << " 2=" << a.Get(256) << " -4=" << a.Get(348) 
       << " 12=" << a.Get(865) << " -25=" << a.Get(1321) << endl;


  a.Set(110,356);
  cerr << "356=" << a.Get(110) << " 64=" << a.Get(18) << " 2=" << a.Get(256) << " -4=" << a.Get(348) 
       << " 12=" << a.Get(865) << " -25=" << a.Get(1321) << endl;
  
  
  cerr << "Now we take up " << a.MemoryUsage() << " bytes" << endl;
  cerr << "done" << endl;


  Array b(100000,0);

  cerr << "sparse array of " << b.size() << " elements;" << endl;
  cerr << "takes up " << b.MemoryUsage() << " bytes;\nStarting indexing operators test:"<< endl;

  b[2]=2;
  b[4]=4;
  b[8]=8;
  b[9]=9;
  b[256]=2;
  b[348]=-4;
  b[865]=12;
  b[1321]=-25;
  b[2]+=3;
  cerr << "5=" << b[2] << " 4=" << b[4] 
       << " 8=" << b[8] << " 9=" << b[9] << endl;
  
  cerr << "2=" << b[256] << " -4=" << b[348] 
       << " 12=" << b[865] << " -25=" << b[1321] << endl;

  cerr << "Breaking the sortedness: " << endl;

  b[18]=64;
  cerr << "64=" << b[18] << " 2=" << b[256] << " -4=" << b[348] 
       << " 12=" << b[865] << " -25=" << b[1321] << endl;


  b[110]=356;
  cerr << "356=" << b[110] << " 64=" << b[18] << " 2=" << b[256] << " -4=" << b[348] 
       << " 12=" << b[865] << " -25=" << b[1321] << endl;
  

  cerr << "Now we take up " << b.MemoryUsage() << " bytes" << endl;
  cerr << "done" << endl;

  
  cerr << "Timing access:" << endl;

#define TESTARRAYSIZE 1000000
  Array c(TESTARRAYSIZE,0);
  c.ReserveStorage(20000);
  vector<int> v(20000);
  TaskTimer t;
  
  cerr << "vector of int, pre-reserved, 1000 cycles, setting  20000 elements: " << endl;
  t.Start();
  for ( int cycle = 0 ; cycle < 1000 ; cycle++ ) {
    for ( unsigned int i = 0 ; i < 20000 ; i++ ) {
      v[i] = i+1;
    }
  }
  t.Stop();
  cerr << "set value timing: " << t << endl;

  t.Reset();
  cerr << "Sorted sparse vector; 1000 cycles, setting 20000 elements by indexing without reallocation: " << endl;
  t.Start();
  for ( int cycle = 0 ; cycle < 1000 ; cycle++ ) {
    for ( unsigned int i = 0 ; i < 20000 ; i++ ) {
      c[i] = i+1;
    }
  }
  t.Stop();
  cerr << "set value timing: " << t << endl;

  t.Reset();
  cerr << "Sorted sparse vector; 1000 cycles, setting 20000 elements by setter without reallocation: " << endl;
  t.Start();
  for ( int cycle = 0 ; cycle < 1000 ; cycle++ ) {
    for ( unsigned int i = 0 ; i < 20000 ; i++ ) {
      c.Set(i,i+1);
    }
  }
  t.Stop();
  cerr << "set value timing: " << t << endl;


  t.Reset();
  cerr << "Sorted sparse vector; 1000 cycles, setting 20000 elements by indexing with erasing/reallocation: " << endl;
  Array d(TESTARRAYSIZE,0);
  t.Start();
  for ( int cycle = 0 ; cycle < 1000 ; cycle++ ) {
    d.ResetToDefault();
    for ( unsigned int i = 0 ; i < 20000 ; i++ ) {
      d[i] = i+1;
    }
  }
  t.Stop();
  cerr << "set value timing: " << t << endl;


  t.Reset();
  cerr << "Sorted sparse vector; 1000 cycles, setting 20000 elements by setter with erasing/reallocation: " << endl;
  t.Start();
  for ( int cycle = 0 ; cycle < 1000 ; cycle++ ) {
    d.ResetToDefault();
    for ( unsigned int i = 0 ; i < 20000 ; i++ ) {
      d.Set(i,i+1);
    }
  }
  t.Stop();
  cerr << "set value timing: " << t << endl;

  t.Reset();
  cerr << "Sorted sparse vector; 1000 cycles, retrieving 20000 elements by non-const indexing: " << endl;
  Array e(TESTARRAYSIZE,0);
  for ( unsigned int i = 0 ; i < 20000 ; i++ ) {
    e.Set(i,i);
  }
  t.Start();
  for ( int cycle = 0 ; cycle < 1000 ; cycle++ ) {
    for ( unsigned int i = 0, j=10000 ; i < 20000 ; i++, j++ ) {
      v[i] = e[j];//d.Get(j);
    }
  }
  t.Stop();
  cerr << "set value timing: " << t << endl;


  t.Reset();
  cerr << "Sorted sparse vector; 1000 cycles, retrieving 20000 elements by const indexing: " << endl;
  const Array & f = e;
  t.Start();
  for ( int cycle = 0 ; cycle < 1000 ; cycle++ ) {
    for ( unsigned int i = 0, j=10000 ; i < 20000 ; i++, j++ ) {
      v[i] = f[j];//d.Get(j);
    }
  }
  t.Stop();
  cerr << "set value timing: " << t << endl;


  t.Reset();
  cerr << "Sorted sparse vector; 1000 cycles, retrieving 20000 elements by iterator: " << endl;
  t.Start();
  v.resize(TESTARRAYSIZE);
  for ( int cycle = 0 ; cycle < 1000 ; cycle++ ) {
    Array::iterator it = e.begin();
    for ( int i = 0 ; i < 20000 ;++it, ++i ) {
      v[i] = *it;
    }
  }
  t.Stop();
  cerr << "set value timing: " << t << endl;
  
  cerr << "Sparse array of 20,000 elements takes up " << e.MemoryUsage()/1000.0 << " Kbytes" << endl;

  cerr << "Checking iterator increment/decrement operators:\n";

  Array::const_iterator iter = a.begin();
  int i = 0;
  for ( ; i <= 10 ; i++ ) {
    cerr << *iter << " ";
    iter++;
  }
  cerr << "\nRolling back to position 2 and starting over:\n" ;
  for ( ; i > 2 ; i-- ) {
    --iter;
  }
  for ( int j = 0 ; j < i ; j++ ) cerr << "  ";
  for ( ; i <= 10 ; i++ ) {
    cerr << *iter << " ";
    iter++;
  }
  cerr << endl;
}
