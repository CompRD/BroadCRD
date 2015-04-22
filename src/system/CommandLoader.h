// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

// CommandLoader manages subprocesses.

#ifndef COMMANDLOADER_H
#define COMMANDLOADER_H

#include "CoreTools.h"

#include <sys/wait.h>
#include <set>

using std::set;

typedef int CommandFunction(int, char**);

class CommandLoaderSingletonGuard;
class CommandLoader
{
    friend class CommandLoaderSingletonGuard;

  public:
    class WaitInfo
    {
        // let CommandLoader have direct access into WaitInfo
        friend class CommandLoader;

      public:
        WaitInfo() : m_pid(-1), m_status(-1), m_timeOut(-1) {}
        WaitInfo(pid_t pid, int waitTime = -1) : m_pid(pid), m_status(-1),
                                                 m_timeOut(waitTime) {}

        pid_t GetPID() const { return m_pid; }
        int GetExitCode() const { return WEXITSTATUS(m_status); }
        int GetSignal() const { return WTERMSIG(m_status); }
        bool IsSignaled(int& nSignal)
        {
            if (WIFSIGNALED(m_status))
            {
                nSignal = WTERMSIG(m_status);
                return true;
            }
            return false;
        }

      private:
        pid_t m_pid;
        int   m_status;
        int   m_timeOut;
    };

  public:
    virtual ~CommandLoader() { StopAllRuns(); }

    // Run() runs the command in strCommand by forking a new AssembleTimer
    // process, which will then in turn launches the command.  The optional
    // parameter will redirect the output of command to either a file or
    // a file descriptor.
    // Run() returns the pid of the AssembleTimer process if successful.
    // Otherwise it returns -1
    // If append=true, RUN appends to the log file
    virtual pid_t Run(const String& strCommand, 
                      const String& strLogFile,
                      const bool useTimer = true,
                      const bool appendToLog = false);

    virtual pid_t Run(const String& strCommand, 
                      int fileDescriptor,
                      const bool useTimer = true );

    virtual pid_t Run(const String& strCommand)
    { 
        return Run(strCommand, STDOUT_FILENO); 
    }

    // the following set of Run() runs a function instead of an executable
    virtual pid_t Run(CommandFunction* pFuncPtr, const String& strCommand,
                      const String& strLogFile = "/dev/null",
                      bool bAppend = false);
    virtual pid_t Run(CommandFunction* pFuncPtr, const String& strCommand,
                      int fdOut);

    // Wait() waits for a child process to terminate.  If a child pid is
    // specified, it will wait for that child to complete.
    // Return values:
    // pid >  0: pid of the child process that exited
    // pid == 0: no child process has exited.  Only valid with a timeout value
    //               of 0.
    // pid = -1: an unknown error has occurred.
    pid_t Wait(int* nExitStatus = NULL, pid_t pid = -1, int timeout = -1);

    // Wait() waits for a child process to terminate.  Optionally, a PID can
    // be specified within the WaitInfo structure, which will cause Wait()
    // to wait for the specified PID.  If successful, Wait() returns true,
    // and the waitInfo structure will be filled with information regarding
    // the exited process.  If Wait() encountered an error, it will return
    // false and waitInfo's data is undefined
    bool  Wait(WaitInfo& waitInfo);

    // singleton access functions
    static CommandLoader* GetInstance();
    // destroy the singleton
    static void Destroy();
    // destroy the singleton and exit with the specified exit code
    static void DestroyAndExit(int nExitCode) { Destroy(); exit(nExitCode); }

    // stops all currently running runs.  StopAllRuns() will not return
    // until all runs have exited
    void StopAllRuns(int signo = SIGINT);

    // send the specified signal to all runs and return immediately
    void SignalAllRuns(int signo) const;

    // return true if CommandLoader is waiting for the exit status of
    // processes
    bool IsRunning() const { return m_setPIDCommands.size() > 0; }
    static bool DisplayExitMessage() { return s_bDisplayExitMessage; }

    virtual void Log(const String& strMessage) const;

  protected:
    // protected constructor to enforce singleton instance
    CommandLoader();
    static CommandLoader* s_pLoader;
    static bool           s_bDisplayExitMessage;

  private:
    void InstallSignalHandlers();
    void RestoreSignalHandlers();

    String m_strAssembleTimerPath;

    set<pid_t> m_setPIDCommands;
    struct sigaction m_oldSigIntAction;
    struct sigaction m_oldSigTermAction;
};

class AssembleCommandLoader : public CommandLoader
{
  public:
    static void Initialize(const String& strLogFile)
    {
        delete s_pLoader;
        s_pLoader = new AssembleCommandLoader(strLogFile);
        s_bDisplayExitMessage = true;
    }

    static void Destroy()
    {
        delete s_pLoader;
        s_pLoader = NULL;
    }

    ~AssembleCommandLoader() { ShutDown(); }

    // Because one of the virtual functions of Run() has been overloaded here,
    // all five of them have to be overloaded, to avoid a compiler warning.

    virtual pid_t Run(const String& strCommand, 
		      const String& strLogFile,
		      const bool useTimer = true,
		      const bool appendToLog = false)
    { return CommandLoader::Run(strCommand, strLogFile, useTimer, appendToLog); }

    virtual pid_t Run(const String& strCommand, 
                      int fileDescriptor,
                      const bool useTimer = true )
    { return CommandLoader::Run(strCommand, fileDescriptor, useTimer); }

    // This is the only Run() method that behaves differently than
    // that provided by CommandLoader.
    virtual pid_t Run(const String& strCommand)
    { return CommandLoader::Run(strCommand,  m_fdWrite); }

    virtual pid_t Run(CommandFunction* pFuncPtr,
		      const String& strCommand,
                      const String& strLogFile = "/dev/null",
                      bool bAppend = false)
    { return CommandLoader::Run(pFuncPtr, strCommand, strLogFile, bAppend); }

    virtual pid_t Run(CommandFunction* pFuncPtr,
		      const String& strCommand,
                      int fdOut)
    { return CommandLoader::Run(pFuncPtr, strCommand, fdOut); }


  bool ShutDown();
    void Log( const String &message ) const;

  private:
    // private constructor to ensure that CommandLoader is a singleton.
    AssembleCommandLoader(const String& strLogFile);

    // create a child process to handle the contents of the message pipe by
    // printing it to stdout and writing it out to a log file
    pid_t LaunchTail(const String& strLogFile);

  private:
    pid_t   m_pidTail;
    int     m_fdRead;
    int     m_fdWrite;
};

#endif
