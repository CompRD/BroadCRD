///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * TestParallelForEach.cc
 *
 *  Created on: Jul 31, 2014
 *      Author: tsharpe
 */

#include "IteratorRange.h"
#include "feudal/ParallelForEach.h"
#include "system/SysConf.h"
#include <atomic>
#include <iostream>

int main( int argc, char** argv )
{
    //configNumThreads(1);
    std::atomic_size_t sum(0);
    parallelForEach(rangeItr(1ul),rangeItr(1001ul),
                    [&sum]( size_t val ){ sum += val; });
    std::cout << "The sum of the first thousand integers is " << sum << std::endl;
}
