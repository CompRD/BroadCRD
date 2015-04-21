///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "system/Pthread.h"
#include <unistd.h>

void Thread::caller( void* ppp )
{
    Thread* pThread = static_cast<Thread*>(ppp);
    pThread->_state = RUNNING;
    pThread->_returnvalue = (*pThread->_threadfunc)(pThread->_argument);
    pThread->_state = COMPLETE;
}

void Thread::Execute( void* arg )
{
    if ( arg )
        _argument = arg;

    if ( !_threadfunc )
        FatalErr("Thread error: no function to execute was ever bound to the "
                 "thread.  Call the Thread::Bind(...) method first.");

    _thread = std::thread(caller,this);
}

void* Thread::Wait()
{
    if ( _thread.joinable() )
        _thread.join();
    return _returnvalue;
}

void vecThread::Run( unsigned jobs )
{
    unsigned int totaljobs = size();
    unsigned int completedjobs = 0;
    unsigned int concurrentjobs = (!jobs || jobs>totaljobs) ? totaljobs : jobs;
    unsigned int currentjob = 0;
    unsigned int runningjobs = 0;

    while ( completedjobs < totaljobs )
    {
        while ( runningjobs < concurrentjobs && currentjob < totaljobs )
        {
            (*this)[currentjob].Execute();
            currentjob++;
            runningjobs++;
        }

        do
        {
            sleep(1);
            for ( Thread& thread : *this )
            {
                if ( thread._state == Thread::COMPLETE )
                {
                    thread.Wait();
                    thread._state = Thread::CLEARED;
                    completedjobs++;
                    runningjobs--;
                }
            }
        }
        while ( runningjobs == concurrentjobs );
    }
}
