///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * FosmidifyEvents.cc
 *
 *  Created on: Sep 13, 2013
 *      Author: tsharpe
 */
#include "paths/long/varcomp/FosmidifyEvents.h"
#include <algorithm>


// moved as a template to the header file

/*

namespace
{

typedef vec<VariantCall>::iterator VCItr;
bool distinctLoci( VCItr itr, VCItr end )
{
    if ( itr==end ) return true;
    size_t prevOff = itr->getRefOffset();
    while ( ++itr != end )
    {
        if ( itr->getRefOffset()==prevOff )
            return false;
        prevOff = itr->getRefOffset();
    }
    return true;
}

}
void fosmidifyEvents( bvec const& refSeq, vec<VariantCall> const& fosmidEvents,
                        vec<VariantCall>& callerEvents, size_t maxSep,
                        bool verbose )
{
    // make sure no reallocations are necessary so that iterators are not invalidated
    callerEvents.reserve(callerEvents.size()+fosmidEvents.size());

    auto fItr1 = fosmidEvents.begin();
    auto fEnd = fosmidEvents.end();
    auto cItr1 = callerEvents.begin();
    auto cEnd = callerEvents.end();
    while ( fItr1 != fEnd && cItr1 != cEnd )
    {
        size_t prevOffset = std::min(fItr1->getRefOffset(),
                                     cItr1->getRefOffset());
        auto fItr2 = fItr1;
        auto cItr2 = cItr1;
        bool more = true;
        while ( more )
        {
            more = false;
            while ( fItr2!=fEnd && fItr2->getRefOffset() <= prevOffset+maxSep )
            {
                prevOffset = std::max(prevOffset,fItr2->getRefOffset());
                ++fItr2;
                more = true;
            }

            while ( cItr2!=cEnd && cItr2->getRefOffset() <= prevOffset+maxSep )
            {
                prevOffset = std::max(prevOffset,cItr2->getRefOffset());
                ++cItr2;
                more = true;
            }
        }

        // at this point fItr1..fItr2 and cItr1..cItr2 delimit 2 sets of possibly
        // equivalent events within a region.  so check to see if they really
        // are equivalent.
        if ( fItr1 != fItr2 && cItr1 != cItr2 && distinctLoci(cItr1,cItr2) &&
                (fItr2-fItr1 != cItr2-cItr1 || !std::equal(cItr1,cItr2,fItr1)) )
        {
            size_t refStart = std::min(fItr1->getRefOffset(),
                                        cItr1->getRefOffset());
            size_t refStop = std::max(fItr2[-1].getRefOffset()+fItr2[-1].getRefSeqLen(),
                                     cItr2[-1].getRefOffset()+cItr2[-1].getRefSeqLen());
            if ( verbose )
                std::cout << "Testing region " << refStart << " to " << refStop << std::endl;

            bvec fosmidSeq(refSeq.begin()+refStart,refSeq.begin()+refStop);
            bvec callerSeq(fosmidSeq);

            auto fItr3 = fItr2;
            while ( fItr3 != fItr1 )
            {
                --fItr3;
                auto refEvent = fosmidSeq.begin()+(fItr3->getRefOffset()-refStart);
                fosmidSeq.erase(refEvent,refEvent+fItr3->getRefSeqLen());
                bvec const& altSeq = fItr3->getAltSeq();
                bvec newSeq = bvec(fosmidSeq.begin(),refEvent);
                newSeq.append(altSeq);
                newSeq.append(refEvent,fosmidSeq.end());
                fosmidSeq = newSeq;
            }

            auto cItr3 = cItr2;
            while ( cItr3 != cItr1 )
            {
                --cItr3;
                if ( cItr3->isHomozygous() )
                {
                    auto refEvent = callerSeq.begin()+(cItr3->getRefOffset()-refStart);
                    callerSeq.erase(refEvent,refEvent+cItr3->getRefSeqLen());
                    bvec const& altSeq = cItr3->getAltSeq();
                    bvec newSeq = bvec(callerSeq.begin(),refEvent);
                    newSeq.append(altSeq);
                    newSeq.append(refEvent,callerSeq.end());
                    callerSeq = newSeq;
                }
            }

            if ( fosmidSeq == callerSeq )
            {
                if ( verbose )
                {
                    std::cout << "Replace these events:" << std::endl;
                    for ( auto itr=cItr1; itr != cItr2; ++itr )
                        if ( itr->isHomozygous() )
                            std::cout << itr->getRefOffset() << ' '
                                        << itr->getRefSeqLen() << ' '
                                        << itr->getAltSeq() << std::endl;
                    std::cout << "With these:" << std::endl;
                    for ( auto itr=fItr1; itr != fItr2; ++itr )
                        std::cout << itr->getRefOffset() << ' '
                                    << itr->getRefSeqLen() << ' '
                                    << itr->getAltSeq() << std::endl;
                }
                cItr3 = cItr2;
                while ( cItr3 != cItr1 )
                {
                    --cItr3;
                    if ( cItr3->isHomozygous() )
                    {
                        callerEvents.erase(cItr3);
                        --cItr2;
                    }
                }
                callerEvents.insert(cItr2,fItr1,fItr2);
                cItr2 += fItr2-fItr1;
            }
        }

        fItr1 = fItr2;
        cItr1 = cItr2;
    }
    sort(callerEvents.begin(),callerEvents.end());
}
*/
