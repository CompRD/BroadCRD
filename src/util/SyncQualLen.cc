///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file SyncQualLen.cc
 * \author tsharpe
 * \date Feb 7, 2012
 *
 * \brief
 */

#include "MainTools.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "feudal/IncrementalWriter.h"
#include "feudal/VirtualMasterVec.h"
#include <algorithm>
#include <cstddef>

int main( int argc, char** argv )
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_String(FASTB);
    CommandArgument_String(QUALB);
    CommandArgument_String(NEWQUALB);
    EndCommandArguments;

    VirtualMasterVec<bvec> fastb(FASTB.c_str());
    VirtualMasterVec<qvec> qualb(QUALB.c_str());
    IncrementalWriter<qvec> newQualb(NEWQUALB.c_str(),fastb.size());

    typedef VirtualMasterVec<bvec>::const_iterator BItr;
    typedef VirtualMasterVec<qvec>::const_iterator QItr;

    size_t nnn = std::min(fastb.size(),qualb.size());
    BItr end = fastb.begin() + nnn;
    QItr qitr = qualb.begin();
    for ( BItr itr(fastb.begin()); itr != end; ++itr, ++ qitr )
    {
        qvec qv = *qitr;
        qv.resize(itr->size());
        newQualb.add(qv);
    }
    BItr realEnd = fastb.end();
    qvec qv;
    while ( end != realEnd )
    {
        qv.resize(end->size());
        ++end;
    }
    newQualb.close();
}
