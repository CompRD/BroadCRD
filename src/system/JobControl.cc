/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "system/CommandLoader.h"
#include "system/JobControl.h"

// initialize JobControl's static variables.
unsigned int JobControl::s_numCPUs = 1;
unsigned int JobControl::s_numJobs = 0;

pid_t JobControl::Run(const String& strCmd, bool* pbAborted, const String* pStrLogFile )
{
    if (!WaitImp(pbAborted) || (pbAborted && *pbAborted))
        return -1;

    String strLogFile;
    
    if ( pStrLogFile ) {
      strLogFile = *pStrLogFile;
    }
    else {
      strLogFile = temp_file::generateName((m_strTmpDir+"/job_XXXXXX").c_str());
    }

    pid_t pid = CommandLoader::GetInstance()->Run(strCmd, strLogFile);

    if (pid == -1)
    {
        cout << "Unable to launch job: " << strCmd << endl;
        return -1;
    }

    m_mapPID2LogFile.insert(make_pair(pid, strLogFile));
    ++s_numJobs;

    cout << "Currently running " << s_numJobs << " jobs on "
         << s_numCPUs << " CPUs" << endl;

    return pid;
}

// Run() runs the specified job by calling the specified command function.
// The notes for the Run() above applies here as well.
pid_t JobControl::Run(CommandFunction* pFunc, const String& strCmd,
                      bool* pbAborted, const String* pStrLogFile )
{
    if (!WaitImp(pbAborted) || (pbAborted && *pbAborted))
        return -1;

    String strLogFile;
    if ( pStrLogFile ) {
      strLogFile = *pStrLogFile;
    }
    else {
      strLogFile = temp_file::generateName((m_strTmpDir + "/job_XXXXXX").c_str());
    }

    pid_t pid = CommandLoader::GetInstance()->Run(pFunc, strCmd,
                                                  strLogFile);

    if (pid == -1)
    {
        cout << "Unable to launch job: " << strCmd << endl;
        return -1;
    }

    cout << "Log file = " << strLogFile << endl;
    m_mapPID2LogFile.insert(make_pair(pid, strLogFile));
    ++s_numJobs;
    return pid;
}

void JobControl::StopAllJobs() const
{
    CommandLoader::GetInstance()->SignalAllRuns(SIGINT);
}

pid_t JobControl::Wait(pid_t pid, const bool keepLogFile )
{
    if (s_numJobs == 0)
    {
        cout << "There are no more jobs" << endl;
        return -1;
    }

    if (pid != -1)  // user requested a specific PID to wait for
    {
        // see if the job has finished already
        bool bFound = false;
        for (list<pid_t>::iterator it(m_listFinishedJobs.begin());
             it != m_listFinishedJobs.end(); ++it)
        {
            if (*it == pid)
            {
                bFound = true;
                m_listFinishedJobs.erase(it);
                break;
            }
        }

        // the job is not finished yet
        if (!bFound)
        {
            CommandLoader::WaitInfo waitInfo(pid);

            if (!CommandLoader::GetInstance()->Wait(waitInfo))
            {
                cout << "Error while waiting for a specific job to finish."
                     << endl;
                return -1;
            }
            else if (waitInfo.GetPID() != pid)
            {
                cout << "I was waiting for PID = " << pid
                     << " but received " << waitInfo.GetPID() << "." << endl;
            }
            pid = waitInfo.GetPID();
            --s_numJobs;
        }
    }
    else if (m_listFinishedJobs.size() > 0)
    {
        pid = m_listFinishedJobs.back();
        m_listFinishedJobs.pop_back();
    }
    else
    {
        pid = CommandLoader::GetInstance()->Wait();

        if (pid == -1)
        {
            cout << "Error while waiting for a job to finish." << endl;
            return -1;
        }
        --s_numJobs;
    }

    map<pid_t, String>::iterator it = m_mapPID2LogFile.find(pid);
    if (it == m_mapPID2LogFile.end())
    {
        cout << "I can't find the log file for the completed job." << endl;
        return -1;
    }

    if (m_pCallback)
    {
        (*m_pCallback)(it->second, m_pData);
    }

    if ( ! keepLogFile )
      unlink(it->second.c_str());

    m_mapPID2LogFile.erase(it);
    return pid;
}

bool JobControl::WaitForAllJobs( const bool keepLogFiles )
{
    while (s_numJobs > 0)
    {
        if (!Wait( -1, keepLogFiles ))
            return false;
    }
    return true;
}

void JobControl::SetNumCPUs(unsigned int numCPUs)
{
    if (numCPUs == 0)
    {
        cout << "Shirley you're not serious when you specified 0 CPUs.  "
             << "Let's pretend you really meant 1." << endl;
        numCPUs = 1;
    }

    s_numCPUs = numCPUs;
}

bool JobControl::WaitImp(bool* pbAborted)
{
    bool bResult = true;

    // All the allotted CPUs are busy.  We'll wait until one of them
    // frees up before continuing
    if (s_numCPUs <= s_numJobs)
    {
        cout << "Waiting for an idle CPU..." << flush;

        CommandLoader::WaitInfo waitInfo;

        if (!CommandLoader::GetInstance()->Wait(waitInfo))
        {
            cout << "Failed to wait for a job to complete" << endl;
            return false;
        }
        else if (pbAborted != NULL)
        {
            int nSignal = 0;
            *pbAborted = (waitInfo.GetExitCode() != 0
                          || waitInfo.IsSignaled(nSignal));
        }

        --s_numJobs;
            // keep the PID of the finished job so that the caller could
        // carry out accounting
        m_listFinishedJobs.push_back(waitInfo.GetPID());

        cout << endl;
    }

    return bResult;
}

// Various convenient callback functions.

void JCCallback_ConcatLogFile(const String& strTmpLogFile, const void* pStrLogFile)
{
    const String& strLogFile = *(static_cast<const String*>(pStrLogFile));
    CpAppend(strTmpLogFile, strLogFile);
}

void JCCallback_PrintLogFile(const String& strTmpLogFile, const void* pData)
{
    System( "cat " + strTmpLogFile );
}

void JCCallback_PrintAndConcatLogFile(const String& strTmpLogFile, const void* pStrLogFile)
{
    JCCallback_ConcatLogFile( strTmpLogFile, pStrLogFile );
    JCCallback_PrintLogFile( strTmpLogFile, 0 );
}

