///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/** \file Pthread.h
This file provides the classes Thread and vecThread - objects representing a
single thread and a thread pool respectively.

The basic steps for spawning a thread are:
 1. Declare an instance of the Thread object
 2. Bind a function to the object
 3. Bind a argument
 4. Execute the Thread
 5. Grab the return value

(For convenience, there are several ways to accomplish these steps, allowing you
 to choose whatever is most appropriate for your program.)

The class vecThread implements a very simple thread pool.  This object stores
 multiple Thread objects and can execute all the threads simultaneously or with
 a user-specified maximum level of concurrency.  In the latter case, a new
 thread starts as soon as an old one finishes (like gmake).  Execution of the
 master thread is blocked until all the spawned threads are complete.

See the file testing/TestPthread.cc for a simple example.
*/

#ifndef _ARACHNE_PTHREAD_H
#define _ARACHNE_PTHREAD_H

#include "Vec.h"
#include <thread>

// Threading function's signature.
typedef void* (*pthrfunc)(void*);

/// Thread: a simple, minimal interface to the POSIX thread library
class Thread
{
public:
    typedef enum { PENDING, RUNNING, COMPLETE, CLEARED } exec_state;

    Thread( pthrfunc fn = nullptr, void* arg = nullptr )
     : _threadfunc(fn), _argument(arg), _returnvalue(nullptr),
       _state(PENDING)
    {}

    /// Binds a function to a thread. Optionally, the argument that should
    /// be passed to the thread can be specified here.
    void Bind( pthrfunc fn, void* arg = nullptr )
    { _threadfunc = fn; _argument = arg; }

    /// Execute a thread.  If an argument has already been specified to the
    /// Bind method, it need not be specified here.
    void Execute( void* arg = nullptr );

    /// Wait until the thread has finished.
    void* Wait();

    /// Returns the pointer to the return value from the threaded function
    void* GetReturnValue() const
    { return _returnvalue; }

private:
    friend class vecThread;
    static void caller( void* );

    pthrfunc _threadfunc;
    void* _argument;
    void* _returnvalue;
    exec_state _state;
    std::thread _thread;
};

/// vecThread: a simple thread pool
class vecThread : public vec<Thread>
{
public:
    vecThread( unsigned nThreads, pthrfunc fn = nullptr )
     : vec<Thread>(nThreads)
    { for ( Thread& thread : *this ) thread.Bind(fn); }

    /// Wait for all of the threads to complete. Useful if you've managed
    /// the thread execution yourself rather than calling Run().
    void Wait()
    { for ( Thread& thread : *this ) thread.Wait(); }

    /// Executes all of the threads in the vecThread object.  Does not wait
    /// for the threads to return (you have to call Wait() explicitly).
    void Execute()
    { for ( Thread& thread : *this ) thread.Execute(); }

    /// Executes all of the threads in the vecThread object.  If jobs is 0,
    /// executes all threads simultaneously.  Otherwise, jobs is used as the
    /// maximum level of concurrency.
    void Run( unsigned jobs = 0 );
};

#endif
