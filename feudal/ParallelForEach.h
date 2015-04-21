///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * ParallelForEach.h
 *
 *  Created on: Jul 31, 2014
 *      Author: tsharpe
 */

#ifndef FEUDAL_PARALLELFOREACH_H_
#define FEUDAL_PARALLELFOREACH_H_

#include "system/Thread.h"
#include "system/SpinLockedData.h"

class ITask
{
public:
    virtual void execute()=0;
    static void parallelExecute( ITask& task );
};

// use a team of threads to process a sequence of values with a functor.
template <class Itr,class Fn>
void parallelForEach( Itr const& beg, Itr const& end, Fn const& fnc )
{
    class Task : public ITask, public SpinLockedData
    {
    public:
        Task( Itr const& beg, Itr const& end, Fn const& fn )
        : mBeg(beg), mEnd(end), mFn(fn) {}

        void execute() override
        { Itr item; Fn fn(mFn);
          while ( next(item) ) fn(*item); }

    private:
        bool next( Itr& item )
        { SpinLocker lock(*this);
          if ( mBeg != mEnd ) { item = mBeg; ++mBeg; return true; }
          return false; }

        Itr mBeg;
        Itr mEnd;
        Fn const& mFn;
    };

    if ( !inParallelSection() )
    {
        Task task(beg,end,fnc);
        ITask::parallelExecute(task);
    }
    else
    {
        Fn fn(fnc);
        for ( Itr itr = beg; itr != end; ++itr )
            fn(*itr);
    }
}


#endif /* FEUDAL_PARALLELFOREACH_H_ */
