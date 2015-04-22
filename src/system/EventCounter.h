////////////////////////////////////////////////////////////////////////////
//                  SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//      This software and its documentation are copyright (2009) by the   //
//  Broad Institute.  All rights are reserved.  This software is supplied //
//  without any warranty or guaranteed support whatsoever. The Broad      //
//  Institute is not responsible for its use, misuse, or functionality.   //
////////////////////////////////////////////////////////////////////////////
/*
 * \file EventCounter.h
 * \author tsharpe
 * \date Nov 20, 2009
 *
 * \brief A class to count multi-threaded events.
 * You can use this, for example, to wait until all the items in a Worklist have
 * completed:  The thread procedure can increment the EventCounter when it
 * completes a work-item, and the main thread can wait until the counter reaches
 * the total number of items that it stuffed into the worklist.
 */
#ifndef SYSTEM_EVENTCOUNTER_H_
#define SYSTEM_EVENTCOUNTER_H_

#include "system/LockedData.h"

/// Multi-threaded event counter.
class EventCounter
{
public:
    EventCounter();
   ~EventCounter();

   /// Inform the counter of an event.  Returns the number of events so far.
   unsigned int incrementEventCount();

   /// Reset the event count.
   unsigned int resetEventCount();

   /// Wait until nEvents events have been counted.
   /// Returns the number of events seen.
   unsigned int wait( unsigned int nEvents );

   /// Wait until nEvents (set in the constructor) have been counted or until
   /// a time-out period has expired.
   /// Returns the number of events seen.
   unsigned int wait( unsigned int nEvents, long nSecs );

private:
   EventCounter( EventCounter const& ); // unimplemented -- no copying
   EventCounter& operator=( EventCounter const& ); // unimplemented -- no copying

   unsigned int mCurEvents; // the number of events that have occurred
   unsigned int mMaxEvents; // the number of events we're waiting for
   LockedData mLock; // a mutex
   Condition mCond; // a condition variable for the predicate "nEvents have occurred"
};

#endif /* SYSTEM_EVENTCOUNTER_H_ */
