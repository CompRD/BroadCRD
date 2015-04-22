///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * CountPFQ30.cc
 *
 *  Created on: May 27, 2014
 *      Author: tsharpe
 */
#include "MainTools.h"
#include "lookup/SAM.h"
#include <iostream>

int main( int argc, char** argv )
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_Bool_OrDefault_Doc(RAW, False,
            "use raw counts (all reads; not just PF)");
    CommandArgument_String(BAM);
    EndCommandArguments;

    SAM::BAMFile bam(BAM);
    SAM::Record rec;
    size_t nQ = 0;
    size_t nHQ = 0;
    while ( bam.nextRecord(rec) )
    {
        nQ += rec.getQualityScores().size();
        if ( RAW || rec.isPF() )
            for ( unsigned char q : rec.getQualityScores() )
                if ( q >= 30 )
                    nHQ += 1;
    }
    std::cout << nHQ << ' ' << nQ << ' ' << 1.*nHQ/nQ << std::endl;
}
