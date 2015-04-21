/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file WorklistTest.cc
 * \author tsharpe
 * \date Jan 16, 2009
 *
 * \brief
 *
 *
 */
#include "system/Worklist.h"
#include "random/Random.h"
#include <unistd.h>
#include <iostream>
#include <string.h>
#include "system/ThreadsafeIO.h"
#include "system/EventCounter.h"

using std::endl;
using std::ostream;

class MyProcessor
{
public:
    MyProcessor( EventCounter& ec ) : mID(pthread_self()), mEventCounter(ec)
    {}

    MyProcessor( MyProcessor const& that )
    : mID(pthread_self()), mEventCounter(that.mEventCounter)
    {}

    void operator()( char* arg )
    {
        ThreadsafeOStream cout(std::cout);
        cout << "Thread " << mID << " processed the arg: " << arg << endl;
        if ( !strcmp(arg,"exit") )
            pthread_exit(0); // let's see if the monitor makes a replacement
        sleep((randomx()>>10)%3);
        mEventCounter.incrementEventCount();
    };

private:
    pthread_t mID;
    EventCounter& mEventCounter;
};


int main( int argc, char** argv )
{
    ThreadsafeOStream cout(std::cout);
    EventCounter ec;
    MyProcessor p(ec); // You need an instance of your processor type.

    Worklist<char*,MyProcessor> worklist( p, 4 );  // Make a worklist and four threads to process it.

    // Optionally, create a monitor that will respawn threads that die.
    worklist.threadDeathMonitor(5); // In this case, the monitor will wake up every 5 seconds to do the dead-thread check.

    worklist.add( argv, argv+argc/2 ); // Add some work to the list, and the threads will get busy with it.

    // You can do this to wait for the worklist to drain, but you don't normally need to.
    // Note that no items in the worklist doesn't mean that the threads have completed all their processing.  It just means
    // that there aren't any "leftover" workitems on which work has not commenced.
    // Nonetheless, it can be useful to wait for the worklist to drain in case you need to manage memory better than simply
    // hand all the workitems to the worklist at once.
    worklist.waitForEmpty();

    /// Wait for argc/2 workitems to finish.
    cout << "Waiting..." << endl;
    ec.wait(argc/2);
    cout << "Done." << endl;

    // Every now and then you can do something like this to monitor your threads' progress.
    // You could also put this code into a signal handler, so that you could cause a report to be produced
    // by a running application whenever you wanted it.  (You'll have to put a pointer to your worklist in a global variable
    // so the signal handler will know where to find it.)
    for ( int iii = 0; iii < 2; ++iii )
    {
        Worklist<char*,MyProcessor>::ProgressMonitor pm = worklist.getProgress(iii);
        cout << "Processed " << pm.mWorkitem << " at " << ctime(&pm.mTime);
    }

    worklist.add( argv+argc/2, argv+argc ); // Add some more work.

    // This is completely optional, but it's a comfort to some.
    // It's recommended that you just let the worklist fall out of scope, rather than calling this method.
    worklist.waitForDone();

    // The worklist will be destroyed here.
    // The worklist destructor waits until all workitems have been completely processed before returning.
    // So it's safe to just let it fall out of scope.
}
