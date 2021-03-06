///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * FosmidifyEvents.h
 *
 *  Created on: Sep 13, 2013
 *      Author: tsharpe
 */

#ifndef PATHS_LONG_VARCOMP_FOSMIDIFYEVENTS_H_
#define PATHS_LONG_VARCOMP_FOSMIDIFYEVENTS_H_

#include "Basevector.h"
#include "Vec.h"

class VariantCall
{
public:
    VariantCall( size_t refOffset, size_t refSeqLen, bvec const& altSeq,
                    bool isHomozygous )
    : mRefOffset(refOffset), mRefSeqLen(refSeqLen), mAltSeq(altSeq),
      mIsHomozygous(isHomozygous)
    {}

    size_t getRefOffset() const { return mRefOffset; }
    size_t getRefSeqLen() const { return mRefSeqLen; }
    bvec const& getAltSeq() const { return mAltSeq; }
    bool isHomozygous() const { return mIsHomozygous; }

    friend bool operator<( VariantCall const& vc1, VariantCall const& vc2 )
    { return vc1.mRefOffset < vc2.mRefOffset; }

    friend bool operator==( VariantCall const& vc1, VariantCall const& vc2 )
    { return vc1.mRefOffset == vc2.mRefOffset &&
            vc1.mRefSeqLen == vc2.mRefSeqLen &&
            vc1.mAltSeq == vc2.mAltSeq &&
            vc1.mIsHomozygous == vc2.mIsHomozygous; }

    friend bool operator!=( VariantCall const& vc1, VariantCall const& vc2 )
    { return !(vc1==vc2); }

private:
    size_t mRefOffset;
    size_t mRefSeqLen;
    bvec mAltSeq;
    bool mIsHomozygous;
};

// Try to figure out small sets of nearby events made by some caller that can
//    be replaced by an equivalent set of events that we made by aligning fosmid
//    to reference.
// refSeq is the complete sequence of the chromosome with respect to which all
//    offsets are given.
// The vector of fosmid events is expected to be sorted, and monotonically
//    increasing on reference offset.
// The vector of caller events is modified in place to make it more like the
//    fosmid calls iff the replacement is completely equivalent to the original
//    calls.
// The vector of caller events is also expected to be sorted, but it may have
//    multiple events at a given reference offset (in which case it won't be
//    considered for fosmidification in that region).

// note (blau): events are either deleted from callerEvents or copied from fosmidEvents

template <typename var_t>
void fosmidifyEvents( bvec const& refSeq, vec<var_t> const& fosmidEvents,
                        vec<var_t>& callerEvents, size_t maxSep=100,
                        bool verbose =false)
{
    typedef typename vec<var_t>::iterator VCItr;
    auto distinctLoci = []( VCItr itr, VCItr end )
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
                          };
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
                (fItr2-fItr1 != cItr2-cItr1 || !std::equal(cItr1,cItr2,fItr1
                                                          ,[](const var_t&L,const var_t&R)
                                                             { return    L.getRefOffset()==R.getRefOffset()
                                                                      && L.getRefSeqLen()==R.getRefSeqLen()
                                                                      && L.getAltSeq()==R.getAltSeq()
                                                                      && L.isHomozygous()==R.isHomozygous();
                                                             }
                                                          )
                )
           )
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

#endif /* PATHS_LONG_VARCOMP_FOSMIDIFYEVENTS_H_ */
