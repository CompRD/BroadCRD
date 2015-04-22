///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// ForkManager Class.
// Dispatches jobs using Fork, in parallel if requested.
// 
// 1) Create a ForkManager object.
// 2) Build a queue using the add method
// 3) Run the jobs using the run method
// 
//    ForkManager myqueue;
//    myqueue.add("echo job1");
//    myqueue.add("echo job2");
//    int failed = myqueue.run();
//
// Run multiple jobs at once by setting the concurrency level > 1
//   ForkManager myqueue(3);


#include "String.h"
#include "Vec.h"

#ifndef PAIRS_MANAGER_H
#define PAIRS_MANAGER_H

class ForkManager {

public:
  enum State {QUEUED, RUNNING, FINISHED, FAILED};

private:
  const vec<String> state_str {"QUEUED", "RUNNING", "FINISHED", "FAILED"};
  uint32_t m_concurrency;
  vec<String> m_commands;
  vec<State> m_job_status;
  
public:
  ForkManager() : m_concurrency(1) {};
  ForkManager(uint32_t concurrency) : m_concurrency(concurrency) {};
  
  // Add job to the queue
  void add(const String& command) ;
  
  // Run all jobs in the queue - returns count of failed jobs
  int run ();
  
  // Reset job status to QUEUED
  void reset() {
    m_job_status.assign(m_job_status.size(), QUEUED);
  }

  // Get a list of job commands
  const vec<String>& get_commands() {
    return m_commands;
  }
  
  // Get a list of job statuses
  const vec<State>& get_status() {
    return m_job_status;
  }
  
  // Overloaded output stream
  friend ostream& operator <<(ostream &os, const ForkManager &obj);
  
};

#endif

