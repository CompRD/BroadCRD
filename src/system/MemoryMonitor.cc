///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * MemoryMonitor.cc
 *
 *  Created on: Mar 22, 2013
 *      Author: tsharpe
 */

#include "system/MemoryMonitor.h"
#include "system/SysConf.h"
#include "system/System.h"
#include <sys/time.h>
#include <sys/resource.h>
#include <ctime>
#include <utility>

MemoryMonitor::MemoryMonitor( unsigned pollSecs, char const* logFile )
  : mPollSecs(pollSecs), mStatm("/proc/self/statm"), mTimeBase(time(nullptr)),
    mCurrentRSS(0), mBaseRSS(0), mMaxRSS(0), mDone(false)
{
    if ( logFile )
    {
        mLog.open(logFile,std::ios_base::out|std::ios_base::binary|std::ios_base::trunc);
        if ( !logFile )
            FatalErr("Unable to open memory monitor log file: " << logFile);
        mLog.write(reinterpret_cast<char*>(&mTimeBase),sizeof(mTimeBase));
        startThread();
    }
}

MemoryMonitor::~MemoryMonitor()
{
    if ( mThread.joinable() )
    {
        mDone = true;
        mCondVar.notify_one();
        mThread.join();
    }
    if ( mLog.is_open() )
    {
        setMaxRSS();
        unsigned maxRSSMb = mMaxRSS/(1024ul*1024ul);
        std::pair<unsigned,unsigned> rec(0u,maxRSSMb);
        mLog.write(reinterpret_cast<char*>(&rec),sizeof(rec));
        mLog.close();
    }
}

void MemoryMonitor::setCurRSS()
{
    mStatm.seekg(0);
    unsigned vmsize, vmrss;
    mStatm >> vmsize >> vmrss;
    mCurrentRSS = vmrss * pageSize();
}

void MemoryMonitor::setMaxRSS()
{
    struct rusage ruse;
    getrusage(RUSAGE_SELF,&ruse);
    mMaxRSS = ruse.ru_maxrss * 1024ul;
}

void MemoryMonitor::pollRSS()
{
    while ( !mDone )
    {
        std::unique_lock<std::mutex> lock(mMutex);
        if ( mCondVar.wait_for(lock,std::chrono::seconds(mPollSecs)) ==
                std::cv_status::timeout )
        {
            setCurRSS();
            if ( mCurrentRSS > mMaxRSS )
                mMaxRSS = mCurrentRSS;
            if ( mLog.is_open() )
            {
                unsigned now = time(nullptr) - mTimeBase;
                unsigned curRSSMb = mCurrentRSS/(1024ul*1024ul);
                std::pair<unsigned,unsigned> rec(now,curRSSMb);
                mLog.write(reinterpret_cast<char*>(&rec),sizeof(rec));
            }
        }
    }
}
