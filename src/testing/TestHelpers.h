/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


/// Short, inline helper functions for testing.
/// \file TestHelpers.h

#ifndef TEST_HELPERS_H
#define TEST_HELPERS_H

#include "String.h"
#include "random/Random.h"

inline void StringSetRandom(String & s, int size) {
  if (s.size() != (size_t) size) { s.resize((size_t) size); }
  for (int i = 0; i != size; ++i) {
    s[i] = randomx() % 26 + 'A';
  }
}

#endif //TEST_HELPERS_H
