////////////////////////////////////////////////////////////////////////////
//                  SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//      This software and its documentation are copyright (2009) by the   //
//  Broad Institute.  All rights are reserved.  This software is supplied //
//  without any warranty or guaranteed support whatsoever. The Broad      //
//  Institute is not responsible for its use, misuse, or functionality.   //
////////////////////////////////////////////////////////////////////////////
/*
 * \file EventCounter.cc
 * \author tsharpe
 * \date Nov 20, 2009
 *
 * \brief
 */
#include "system/EventCounter.h"

EventCounter::EventCounter()
: mCurEvents(0), mMaxEvents(0), mCond(mLock)
{
}

EventCounter::~EventCounter()
{
}

unsigned int EventCounter::incrementEventCount()
{
    Locker lock(mLock);
    if ( ++mCurEvents >= mMaxEvents )
        mCond.signal();
    return mCurEvents;
}

unsigned int EventCounter::resetEventCount()
{
    Locker lock(mLock);
    unsigned int result = mCurEvents;
    mCurEvents = 0;
    return result;
}

unsigned int EventCounter::wait( unsigned int nEvents )
{
    Locker lock(mLock);
    mMaxEvents = nEvents;
    if ( mCurEvents < mMaxEvents )
        lock.wait(mCond);
    return mCurEvents;
}


unsigned int EventCounter::wait( unsigned int nEvents, long nSecs )
{
    Locker lock(mLock);
    mMaxEvents = nEvents;
    if ( mCurEvents < mMaxEvents )
        lock.timedWait(mCond,nSecs);

    return mCurEvents;
}
