#ifndef TEST_ASSERT_H
#define TEST_ASSERT_H

#include "CoreTools.h"
#include "system/Assert.h"

#define LocationInfo cout << "line " << __LINE__ \
                          << ", function " << __FUNCTION__ \
                          << ", file " << __FILE__ << endl;

///Fail if X and Y (of any type) are further apart than EPSILON/mean(|X|,|Y|).
#define TestCloserEpsilon(X,Y,EPSILON) if (!TestCloseFunction(X,Y,EPSILON)) { \
cout << "AssertClose: "; PRINT2(X,Y); LocationInfo; CRD::exit(1); }\
 else {}

///Fail if X and Y (of any type) are further apart than .0001/mean(|X|,|Y|).
#define TestClose(X,Y) if (!TestCloseFunction(X,Y)) { \
cout << "AssertClose: "; PRINT2(X,Y); LocationInfo; CRD::exit(1); }\
 else {}

///False if x and y are further apart than epsilon/mean(|X|,|Y|).
///X and Y must support +, -, and automatic cast to double.
template<class X, class Y>
bool TestCloseFunction(const X & x, const Y & y, double epsilon = 0.0001) {
  double a = abs(x) + abs(y);
  if (a != 0) epsilon *= a/2;
  return abs(x-y) < epsilon;
}

inline void Checkpoint(const String &label = String("Test ")) {
  static int iteration = 0;
  cout << label << ++iteration << endl;
}
#endif //TEST_ASSERT_H
