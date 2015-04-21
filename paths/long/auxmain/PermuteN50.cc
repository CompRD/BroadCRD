///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * PermuteN50.cc
 *
 *  Created on: Jul 25, 2014
 *      Author: tsharpe
 */
#include "MainTools.h"
#include "math/Functions.h"
#include <algorithm>
#include <fstream>
#include <vector>

namespace
{

class N50Calculator
{
public:
    N50Calculator()
    : mPending(0) { mVals.reserve(10000); }

    void add( unsigned val )
    { if ( !val ) return;
      mPending += val; }

    void flush()
    { if ( mPending ) { mVals.push_back(mPending); mPending = 0; } }

    unsigned getN50()
    { flush();
      if ( mVals.empty() ) return 0;
      std::sort(mVals.begin(),mVals.end());
      return N50(mVals); }

private:
    vec<unsigned> mVals;
    unsigned mPending;
    std::ofstream mOS;
};

}

int main( int argc, char** argv )
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_String_Doc(N50, "name of N50 log file");
    EndCommandArguments;

    std::vector<std::vector<unsigned>> vals;
    vals.reserve(1000);
    std::ifstream in(N50);
    int val;
    while ( in.read(reinterpret_cast<char*>(&val),sizeof(val)) )
    {
        if ( val < 0 ) vals.resize(vals.size()+1);
        else vals.back().push_back(val);
    }
    int count = 100;
    while ( count-- )
    {
        std::random_shuffle(vals.begin(),vals.end());
        N50Calculator n50calc;
        for ( auto const& vec : vals )
            for ( unsigned val : vec )
                if ( !val ) n50calc.flush();
                else n50calc.add(val);
        std::cout << "N50=" << n50calc.getN50() << std::endl;
    }
}
