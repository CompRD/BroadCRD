///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * ParallelForEach.cc
 *
 *  Created on: Jul 31, 2014
 *      Author: tsharpe
 */
// MakeDepend: library OMP
// somewhat surprisingly, the next line is necessary even though we're using
// std::thread and we're not using pthreads directly
// MakeDepend: library PTHREAD

#include "feudal/ParallelForEach.h"
#include "system/LockedData.h"
#include "system/SysConf.h"
#include <atomic>
#include <thread>
#include <omp.h>

namespace
{

class ThreadPool
{
public:
    ThreadPool( ThreadPool const& )=delete;
    ThreadPool& operator=( ThreadPool const& )=delete;

    ThreadPool()
    : mNThreads(getConfiguredNumThreads()),
      mpThreads(new std::thread[mNThreads]), mpTask(nullptr), mQuit(false),
      mTaskAvailable(mTaskLock), mTaskComplete(mDoneLock)
    { omp_set_num_threads(1);
      std::thread* pThread = mpThreads+mNThreads;
      while ( pThread-- != mpThreads )
        *pThread = std::thread(&ThreadPool::run,this); }

    ~ThreadPool()
    { if ( true )
      { Locker lockTask(mTaskLock);
        mQuit = true; }
      mTaskAvailable.broadcast();
      std::thread* pThread = mpThreads+mNThreads;
      while ( pThread-- != mpThreads )
        pThread->join();
      omp_set_num_threads(getConfiguredNumThreads()); }

    void execute( ITask& task )
    { if ( true )
      { Locker lockTask(mTaskLock);
        mNDone = 0ul;
        mpTask = &task; }
      mTaskAvailable.broadcast();
      if ( true )
      { Locker lockDone(mDoneLock);
        while ( mNDone != mNThreads )
          lockDone.wait(mTaskComplete); } }

    void run()
    { while ( true )
      { if ( true )
        { Locker lockTask(mTaskLock);
          while ( !mpTask && !mQuit )
            lockTask.wait(mTaskAvailable); }
        if ( mQuit ) break;
        mpTask->execute();
        bool broadcast = false;
        if ( true )
        { Locker lockDone(mDoneLock);
          if ( ++mNDone == mNThreads )
          { mpTask = nullptr; broadcast = true; }
          else
          { while ( mNDone != mNThreads )
              lockDone.wait(mTaskComplete); } }
        if ( broadcast )
          mTaskComplete.broadcast(); } }

private:
    unsigned mNThreads;
    std::thread* mpThreads;

    ITask* mpTask;
    bool mQuit;
    LockedData mTaskLock;
    Condition mTaskAvailable;

    unsigned mNDone;
    LockedData mDoneLock;
    Condition mTaskComplete;
};

}

void ITask::parallelExecute( ITask& task )
{
    static ThreadPool gTP;
    gTP.execute(task);
}
