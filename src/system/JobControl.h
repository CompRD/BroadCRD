/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef __JOB_CONTROL_H_INCLUDED__
#define __JOB_CONTROL_H_INCLUDED__

#include <list>
#include <map>
#include "system/CommandLoader.h"

// this is the prototype for the job control callback used by the
// JobControl object.  When a process exits, the JobControl object will
// invoke this callback by passing the log file associated with that
// process, along with the user specified data pointer.
typedef void JobControlCallback(const String& strLogFile, const void* pData);

class JobControl
{
  public:
    JobControl(const String& strTmpDir, JobControlCallback* pCallback = NULL,
               void* pData = NULL)
        : m_strTmpDir(strTmpDir), m_listFinishedJobs(), m_pCallback(pCallback),
          m_pData(pData)
    {
    }

    ~JobControl() { WaitForAllJobs(); }

    // Run() runs the specified job contained in strCmd.
    // If all of the allotted CPUs are busy with jobs, then Run() will
    // wait until a job finishes.
    // If user specifies a pointer to pbAborted, and if Run() has to wait
    // for a free CPU, Run() will check the exit status of the completed
    // job.  And if the job has exited due to an error or a signal, then
    // Run() will return WITHOUT running the job.
    // If user specifies a pointer to pStrLogFile, output of the
    // command will be placed in that filename.  Otherwise, a
    // temporary filename will be used.
    pid_t Run(const String& strCmd, bool* pbAborted = NULL, 
              const String* pStrLogFile = NULL);

    // Run() runs the specified job by calling the specified command function.
    // The notes for the Run() above applies here as well.
    pid_t Run(CommandFunction* pFunc, const String& strCmd,
              bool* pbAborted = NULL, const String* pStrLogFile = NULL);

    void StopAllJobs() const;
    pid_t Wait(pid_t pid = -1, const bool keepLogFile = false );
    bool WaitForAllJobs( const bool keepLogFiles = false );
    static void SetNumCPUs(unsigned int numCPUs);

  private:
    bool WaitImp(bool* pbAborted);

  private:
    // s_numCPUs and s_numJobs are intentional static members so that you can
    // instantiate and call JobControl anywhere and it will still limit the
    // number of jobs
    static unsigned int s_numCPUs;
    static unsigned int s_numJobs;
    String              m_strTmpDir;
    map<pid_t, String>  m_mapPID2LogFile;
    list<pid_t>         m_listFinishedJobs;
    JobControlCallback* m_pCallback;
    void*               m_pData;
};


// A callback to print the log output from each job to stdout after it
// finishes.  Used like so:
// 
// JobControl jc( tmpDir, &JCCallback_PrintLogFile );

void JCCallback_PrintLogFile(const String& strTmpLogFile, const void* pData);

// A callback to concat all the log output from the jobs to a single
// file.  You must pass in the address of a String containing the name
// of the concatenated log file as the "data" parameter of the
// JobControl constructor, like so:
//
// String masterLogFile( "master_log" );
// JobControl jc( tmpDir, &JCCallback_ConcatLogFile, &masterLogFile );

void JCCallback_ConcatLogFile(const String& strTmpLogFile, const void* pStrLogFile);

// A callback to do both of the above.  Used like JCCallbackConcatLogFile().

void JCCallback_PrintAndConcatLogFile(const String& strTmpLogFile, const void* pStrLogFile);


#endif
