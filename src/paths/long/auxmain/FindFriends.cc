///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file Friends.cc
 * \author tsharpe
 * \date Apr 17, 2012
 *
 * \brief
 */
#include "MainTools.h"
#include "Basevector.h"
#include "paths/long/Friends.h"

int main( int argc, char** argv )
{
    String const EMPTY;
    RunTime();
    BeginCommandArguments;
    CommandArgument_String_Doc(READS, "Reads.");
    CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS,0u,"Threads.");
    CommandArgument_Bool_OrDefault(USE_BROKEN_VERSION,True);
    EndCommandArguments;

    configNumThreads(NUM_THREADS);
    vecbvec reads(READS);
    FindFriends(reads,
                 READS.ReplaceExtension(".fastb",".friends"),
                 USE_BROKEN_VERSION);
}
