///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * MemoryMonitor.h
 *
 *  Created on: Mar 22, 2013
 *      Author: tsharpe
 */

#ifndef SYSTEM_MEMORYMONITOR_H_
#define SYSTEM_MEMORYMONITOR_H_

#include <cstddef>
#include <fstream>
#include <condition_variable>
#include <mutex>
#include <thread>

class MemoryMonitor
{
public:
    MemoryMonitor( unsigned pollSecs = 20, char const* logFile = 0 );
    MemoryMonitor( MemoryMonitor const& ) = delete;
    MemoryMonitor& operator=( MemoryMonitor const& ) = delete;
   ~MemoryMonitor();

    // last observed RSS (in bytes)
    size_t getCurrentRSS() const { return mCurrentRSS; }

    void trackUse()
    { startThread();
      mBaseRSS = mMaxRSS = mCurrentRSS; }

    size_t getUse() { return mMaxRSS - mBaseRSS; }

private:
    void setCurRSS();
    void setMaxRSS();

    void startThread()
    { if ( !mThread.joinable() )
      { setCurRSS();
        mThread = std::thread(&MemoryMonitor::pollRSS,this); } }

    void pollRSS();

    unsigned mPollSecs;
    std::ofstream mLog;
    std::ifstream mStatm;
    long mTimeBase;
    size_t mCurrentRSS;
    size_t mBaseRSS;
    size_t mMaxRSS;
    std::mutex mMutex;
    std::condition_variable mCondVar;
    std::thread mThread;
    bool mDone;
};


#endif /* SYSTEM_MEMORYMONITOR_H_ */
