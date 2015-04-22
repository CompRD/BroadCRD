///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file Test.cc
 * \author tsharpe
 * \date Jun 12, 2011
 *
 * \brief
 */
#include "MainTools.h"
#include "system/WorklistUtils.h"
#include <iostream>

int main( int argc, char *argv[] )
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_LongLong(UNITS);
    CommandArgument_UnsignedInt(BYTES_PER_UNIT);
    CommandArgument_UnsignedInt(MIN_BATCH_SIZE);
    CommandArgument_UnsignedInt(NUM_THREADS);
    CommandArgument_LongLong(MAX_MEMORY_GB);
    EndCommandArguments;
    SetMaxMemory(MAX_MEMORY_GB<<30);

    WorklistParameterizer wp(UNITS,BYTES_PER_UNIT,MIN_BATCH_SIZE,NUM_THREADS);
    std::cout << "Doing " << UNITS << " units of work using "
              << wp.getNThreads() << " threads on " << wp.getNBatches()
              << " batches of " << wp.getBatchSize() << '.' << std::endl;
}
