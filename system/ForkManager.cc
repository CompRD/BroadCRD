///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <sys/wait.h>
#include "system/ForkManager.h"

void ForkManager::add(const String& command) {
    m_commands.push_back(command);
    m_job_status.push_back(QUEUED);
}
  
// Run all jobs in the pool
int ForkManager::run () {
    // Check that there is something to do
  if (m_commands.empty()) {
    cout << "Warning: No jobs found." << endl;
    return 0;
  }
  
  uint32_t job_count = m_commands.size();
  cout << "Queue contains " << job_count << " jobs." << endl;
  
  vec<int> pids(job_count, 0);
  // While queued jobs remain
  while (std::count(m_job_status.begin(), m_job_status.end(), FINISHED ) +
	 std::count(m_job_status.begin(), m_job_status.end(), FAILED ) < job_count) {
    // While free slots remain
    while (std::count(m_job_status.begin(), m_job_status.end(), RUNNING) < m_concurrency &&
	   std::find(m_job_status.begin(), m_job_status.end(), QUEUED) != m_job_status.end()) {
      uint32_t next = std::find(m_job_status.begin(), m_job_status.end(), QUEUED) - m_job_status.begin();
      pids[next] = Fork(m_commands[next]);
      m_job_status[next] = RUNNING;
      cout << "Starting job: " << next << endl;
    }
    
    // Wait for job to complete
    int status = -2;
    int pid = wait( &status );
    
    int index = std::find(pids.begin(), pids.end(), pid) - pids.begin();
    ForceAssert( index >= 0 );
    pids[index] = 0;
    
    // Report job status
    if (status == 0) {
      cout << "Finished job: " << index << endl;
      m_job_status[index] = FINISHED;
    }	else if ( status != 0 ) {
      cout << "Failed job: " << index << endl;
      m_job_status[index] = FAILED;
    }
    
  }
  // Return failed job count
  return std::count(m_job_status.begin(), m_job_status.end(), FAILED );
}


ostream& operator <<(ostream &os, const ForkManager &obj) {
  for (uint32_t index  = 0; index < obj.m_commands.size(); ++index) 
    os << obj.state_str[obj.m_job_status[index]] << ": " << obj.m_commands[index] << endl;
  return os;
}

