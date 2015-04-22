#ifndef __PAIRS_MANAGER_TEST_H
#define __PAIRS_MANAGER_TEST_H

// Test suite for PairsManager class.
// Tests:
// Construction
// File write and read  


#include <cxxtest/TestSuite.h>
#include <unistd.h>
#include "PairsManager.h"

class MyTest : public CxxTest::TestSuite
{
public:
  // Test fixture setup code.
  void setUp()
  {}
  
  // Test fixture cleanup code.
  void tearDown()
  {}
  
  PairsManager buildPairsManager(longlong pairs = 4, int libraries = 2) {
    PairsManager pm(pairs * 2 + 100);
    for (longlong i = 0; i < pairs; ++i)
      pm.addPair( i*2, i*2 + 1, 20 * (i % libraries), 10 * (i % libraries));  
    return pm;
  }

  void verifyPairsManager(const PairsManager& pm, size_t pairs = 4, size_t libraries = 2) {
    TS_ASSERT_EQUALS(pm.nReads(), pairs * 2 + 100); 
    TS_ASSERT_EQUALS(pm.nPairs(), pairs);
    TS_ASSERT_EQUALS(pm.nLibraries(), libraries); 
    for (size_t i = 0; i < pairs; ++i)
      verifyPair(pm, i, i*2, i*2 + 1, 20 * (i % libraries), 10 * (i % libraries));  
  }

  void verifyPair(const PairsManager& pm, const int index,
		  const int id1, const int id2,
		  const int sep, const int dev) {
    TS_ASSERT_EQUALS(pm.ID1(index), id1); 
    TS_ASSERT_EQUALS(pm.ID2(index), id2); 
    TS_ASSERT_EQUALS(pm.sep(index), sep); 
    TS_ASSERT_EQUALS(pm.sd(index), dev); 
  }

  void writePairsManager(const PairsManager& pm, const String& filename) {
    pm.Write(filename);
  }

  PairsManager readPairsManager(const String& filename) {
    PairsManager pm;
    pm.Read(filename);
    return pm;
  }

  void equalPairsManager(const PairsManager& pm1, const PairsManager& pm2) {
    TS_ASSERT_EQUALS(pm1.nPairs(), pm2.nPairs());
    TS_ASSERT_EQUALS(pm1.nReads(), pm2.nReads());
    TS_ASSERT_EQUALS(pm1.nLibraries(), pm2.nLibraries());
    for (size_t i = 0; i < pm1.nPairs(); ++i)
      equalPair(pm1, pm2, i, i);
    for (size_t i = 0; i < pm1.nLibraries(); ++i)
      TS_ASSERT_EQUALS(pm1.getLibraryName(i), pm2.getLibraryName(i));
  }

  void equalPair(const PairsManager& pm1, const PairsManager& pm2,
		  const int index1, const int index2) {
    TS_ASSERT_EQUALS(pm1.ID1(index1), pm2.ID1(index2)); 
    TS_ASSERT_EQUALS(pm1.ID2(index1), pm2.ID2(index2)); 
    TS_ASSERT_EQUALS(pm1.sep(index1), pm2.sep(index2)); 
    TS_ASSERT_EQUALS(pm1.sd(index1), pm2.sd(index2)); 
    TS_ASSERT_EQUALS(pm1.libraryID(index1), pm2.libraryID(index2));
    TS_ASSERT_EQUALS(pm1[index1], pm2[index2]);
  }

  void test_construction() {
    PairsManager pm = buildPairsManager();
    verifyPairsManager(pm);
  }

  void test_readwrite() {
    const String filename("PairsManagerTests.pairs");

    cout << Date() << " building" << endl;
    PairsManager pm_out = buildPairsManager(2000,10);

    cout << Date() << " writing" << endl;
    writePairsManager(pm_out, filename );

    TS_ASSERT(IsRegularFile(filename));

    cout << Date() << " reading" << endl;
    PairsManager pm_in = readPairsManager(filename );

    cout << Date() << " done" << endl;

    verifyPairsManager(pm_in, 2000, 10);
    equalPairsManager(pm_in, pm_out);
    unlink(filename.c_str());
  }

};

#endif // __MYTEST_H
