/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/// Read space separated floats from stdin; output n, sum, median, average devation.
/// Ignore any non-numeric text mixed in with the floats.
/// \file Stats2.cc

#include "MainTools.h"
#include "math/Functions.h"

int main( int argc, char *argv[] )
{  
  RunTime();
  
  vec<float> data;
  float f;
  String junk;
  int line=0;
  while (!cin.eof()) {
    cin >> f;
    ++line;
    if (cin.fail()) { 
      //PRINT5(line, f, int(cin.eof()), int(cin.fail()), int(cin.bad()));
      if (cin.eof() || cin.bad()) break;
      else {
	cin.clear(); //ignore bad data
	cin >> junk;
      }
    }
    else data.push_back(f); //save good data
    //DotMod(cout, data.size(), 1000);
  }
  Sort(data);
  float median = Median(data,0);
  double avgdev = StdDev(data, median);
  NormalDistribution dist = SafeMeanStdev(data);
  cout << "N= " << data.size()
       << " sum= " << Sum(data)
       << " median= " << median
       << " avgdev= " << avgdev
       << endl;
  return 0;
}
  
