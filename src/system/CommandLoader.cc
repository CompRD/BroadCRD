// Copyright (c) 2000-2004 The Broad Institute at MIT and Harvard

#include "system/CommandLoader.h"
#include "CoreTools.h"

#include <cerrno>
#include <sys/wait.h>
#include <sys/types.h>


// CommandLoader is a singleton.  This is where it lives.
CommandLoader* CommandLoader::s_pLoader = NULL;
bool CommandLoader::s_bDisplayExitMessage = false;

// A guard class that will destroy the command loader singleton when
// the executable containing the command loader terminates, ensuring that
// any and all child processes are cleaned up.
class CommandLoaderSingletonGuard
{
  public:
    CommandLoaderSingletonGuard() {}
    ~CommandLoaderSingletonGuard() { delete CommandLoader::s_pLoader; }
};
static CommandLoaderSingletonGuard s_pCommandLoaderGuard;

CommandLoader::CommandLoader() : m_setPIDCommands()
{
  if ( IsRegularFile( "./AssembleTimer" ) )
    m_strAssembleTimerPath = "./AssembleTimer";
  else
    m_strAssembleTimerPath = search_path_for_command("AssembleTimer");

  if ( m_strAssembleTimerPath.empty() )
  {
    cout << "I cannot find AssembleTimer and I need it.  Exiting."
         << endl;
    exit(-1);
  }
}

static void SigIntHandler(int signo)
{
    if (CommandLoader::GetInstance()->IsRunning())
        CommandLoader::GetInstance()->SignalAllRuns(SIGINT);
    else
    {
        cout << endl << Date( ) + ".  "
             << "Interrupt received (perhaps a ctrl-c).  Stopping." << endl
             << "\nGenerating a backtrace..." << endl;
        TracebackThisProcess( );
        exit(-1);
    }
}

static void SigTermHandler(int signo)
{
    if (!CommandLoader::GetInstance()->IsRunning())
    {
        cout << endl << Date( ) + ".  Killed.  Stopping." << endl
             << "Assembly terminated." << endl << endl;
    }

    CommandLoader::GetInstance()->StopAllRuns(SIGTERM);
    cout << flush;
    cerr << flush;

    if (CommandLoader::DisplayExitMessage())
        exit(-1);
    else
        exit(-1);
}

void CommandLoader::InstallSignalHandlers()
{
    struct sigaction act;
    act.sa_flags = 0;

    // install SIGINT handler
    act.sa_handler = &SigIntHandler;
    sigemptyset(&act.sa_mask);
    sigaddset(&act.sa_mask, SIGINT);
    if (sigaction(SIGINT, &act, &m_oldSigIntAction) < 0)
    {
        cerr << "Error installing SIGINT handler" << endl;
        exit(1);
    }

    // install SIGTERM handler
    act.sa_handler = &SigTermHandler;
    sigemptyset(&act.sa_mask);
    sigaddset(&act.sa_mask, SIGTERM);
    if (sigaction(SIGTERM, &act, &m_oldSigTermAction) < 0)
    {
        cerr << "Error installing SIGTERM handler" << endl;
        exit(1);
    }
}

void CommandLoader::RestoreSignalHandlers()
{
    if (sigaction(SIGTERM, &m_oldSigTermAction, NULL) < 0)
    {
        cerr << "Error restoring SIGTERM handler" << endl;
        exit(1);
    }
    if (sigaction(SIGINT, &m_oldSigIntAction, NULL) < 0)
    {
        cerr << "Error restoring SIGINT handler" << endl;
        exit(1);
    }
}

pid_t CommandLoader::Run(CommandFunction* pFuncPtr, const String &strCommand,
                         const String& strLogFile, bool bAppendToLog)
{
    // get the file creation mask
    int mask = umask(0);
    umask(mask);

    // open the log file for writing
    int fdLogFile(-1);

    fdLogFile = open(strLogFile.c_str(),
		     O_WRONLY | O_CREAT | ( bAppendToLog ? O_APPEND : O_TRUNC ),
		     ~mask);

    if (fdLogFile == -1)
    {
        cerr << "Unable to open " << strLogFile << endl;
        return false;
    }

    pid_t pid = Run(pFuncPtr, strCommand, fdLogFile);
    close(fdLogFile);
    return pid;
}

class FDStdOutGuard
{
  public:
    FDStdOutGuard(int fdNew)
    {
        m_fdStdOutOld = dup(STDOUT_FILENO);
        m_fdStdErrOld = dup(STDERR_FILENO);
        dup2(fdNew, STDOUT_FILENO);
        dup2(fdNew, STDERR_FILENO);
    }

    ~FDStdOutGuard()
    {
        dup2(m_fdStdOutOld, STDOUT_FILENO);
        dup2(m_fdStdErrOld, STDERR_FILENO);
    }

  private:
    int m_fdStdOutOld;
    int m_fdStdErrOld;
};

pid_t CommandLoader::Run(CommandFunction* pFuncPtr, const String& strCommand,
                         int fdOut)
{
    FDStdOutGuard fdStdOutGuard(fdOut);

    pid_t pid = fork();

    if (pid == -1)
        return -1;

    if (pid == 0)
    {
        // make a copy of the command string so that we can transform it
        // into a vector of parameters for execv()
        void* pszCommandBuf = malloc(strCommand.size() + 1);
        char* pszCommand = static_cast<char*>(pszCommandBuf);
        strcpy(pszCommand, strCommand.c_str());

        const char pszDelimiter[] = " ";
        vector<char*> vecArgs;
        char* pszToken;

        char functionLoader[] = "FunctionLoader";
        vecArgs.push_back(functionLoader);
        while ((pszToken = strtok(pszCommand, pszDelimiter)) != NULL)
        {
            vecArgs.push_back(pszToken);
            pszCommand = NULL;
        }

        void*  argvBuf = malloc(sizeof(char*) * vecArgs.size() + 1);
        char** argv = static_cast<char**>(argvBuf);

        for (unsigned int i = 0; i < vecArgs.size(); ++i)
        {
            argv[i] = vecArgs[i];
        }
        argv[vecArgs.size()] = (char*) 0;

        int nResult = (*pFuncPtr)(vecArgs.size(), argv);

	free(argvBuf);
	free(pszCommandBuf);

        _exit(nResult);
    }

    m_setPIDCommands.insert(pid);

    // prevent duplicate signal handler installations
    if (m_setPIDCommands.size() == 1)
        InstallSignalHandlers();

    return pid;
}

pid_t CommandLoader::Run(const String &strCommand,
                         const String& strLogFile,
                         const bool useTimer,
                         const bool appendToLog )
{
    if (strLogFile.size() == 0)
        return -1;

    // get the file creation mask
    int mask = umask(0);
    umask(mask);

    // open the log file for writing
    int fdLogFile(-1);

    fdLogFile = open(strLogFile.c_str(),
		     O_WRONLY | O_CREAT | ( appendToLog ? O_APPEND : O_TRUNC ),
		     ~mask);

    if (fdLogFile == -1)
    {
        cerr << "Unable to open " << strLogFile << endl;
        return false;
    }

    pid_t pid = Run(strCommand, fdLogFile, useTimer );
    close(fdLogFile);
    return pid;
}

pid_t CommandLoader::Run(const String& strCommand, int fdOut,
                         const bool useTimer )
{
    FDStdOutGuard fdStdOutGuard(fdOut);

    pid_t pid = fork();

    if (pid == -1)
        return -1;

    if (pid == 0)
    {

        String strExec, strArgs;
	if ( strCommand.Contains( " " ) ) {
	  strExec = strCommand.Before( " " );
	  strArgs = " " + strCommand.After( " " );
	} else {
	  strExec = strCommand;
	}

        // If the command contains a slash, we'll presume that the
        // caller has specified the path to the command.  Otherwise,
        // we will search the path for the given command.
        String strFullPathOfExec = ( strExec.Contains("/")
                                     ? strExec
                                     : search_path_for_command( strExec ) );

	if ( strFullPathOfExec.empty() ) {
	  cout << "Unable to find " << strExec << " in $PATH." << endl;
	  String PATH( getenv( "PATH" ) );
	  cout << "$"; PRINT( PATH );
	  TracebackThisProcess();
	}

	String strFullCommand = strFullPathOfExec + strArgs;

        // make a copy of the command string so that we can transform it
        // into a vector of parameters for execv()
        void* pszCommandBuf = malloc(strFullCommand.size() + 1);
        char* pszCommand = static_cast<char*>(pszCommandBuf);
        strcpy(pszCommand, strFullCommand.c_str());

        const char pszDelimiter[] = " ";
        vector<char*> vecArgs;
        char* pszToken;

        if ( useTimer )
            vecArgs.push_back(strdup(m_strAssembleTimerPath.c_str()));

        while ((pszToken = strtok(pszCommand, pszDelimiter)) != NULL)
        {
            vecArgs.push_back(pszToken);
            pszCommand = NULL;
        }

        char** argv = static_cast<char**>(malloc(sizeof(char*) * vecArgs.size() + 1));

        for (unsigned int i = 0; i < vecArgs.size(); ++i)
        {
            argv[i] = vecArgs[i];
        }
        argv[vecArgs.size()] = (char*) 0;

	// Note: there is some memory leak here,
	// since pszCommandBuf is not free()ed.

        execv(argv[0], argv);
        cerr << "Unable to run " << argv[0] << ".  Exiting..." << endl;
        exit(-1);
    }

    m_setPIDCommands.insert(pid);

    // prevent duplicate signal handler installations
    if (m_setPIDCommands.size() == 1)
        InstallSignalHandlers();

    return pid;
}

// a signal handler that doesn't do anything, which in essence ignores any
// signals that it is assigned to handle
static void SignalHandler(int signo)
{
}

pid_t CommandLoader::Wait(int* pExitStatus, pid_t pid, int timeout)
{
    if (m_setPIDCommands.size() == 0)
        return -1;

    if (pid != -1) // user specified a specific pid to wait on
    {
        // make sure that this pid belongs to one of our children
        if (m_setPIDCommands.find(pid) == m_setPIDCommands.end())
            return -1;
    }

    int status;
    pid_t p = pid;
    struct sigaction oldAction;

    // if the caller specified a timeout value, set the alarm to go off in
    // this time.
    if (timeout > 0)
    {
        struct sigaction alarmAction;
        alarmAction.sa_handler = &SignalHandler;
        sigemptyset(&alarmAction.sa_mask);
        alarmAction.sa_flags = 0;

        if (sigaction(SIGALRM, &alarmAction, &oldAction) < 0)
        {
            cerr << "Error catching SIGALRM." << endl;
            exit(-1);
        }
        alarm(timeout);
    }

    do
    {
        pid = waitpid(p, &status, timeout == 0 ? WNOHANG : 0);

        // check if we're interrupted by the alarm
        if (pid == -1 && errno == EINTR && timeout > 0)
        {
            // see how much time is left in the alarm.
            unsigned int unslept = alarm(0);

            // if the value is non-zero, that means we were interrupted for
            // some other reason.  Go back to waiting
            if (unslept > 0)
                alarm(unslept);
            else  // the alarm went off, time to exit and return to the caller
                break;
        }
    } while (pid == -1 && errno == EINTR);

    // turn off the alarm
    if (timeout > 0)
    {
        alarm(0);

        // restore the old signal handler
        if (sigaction(SIGALRM, &oldAction, NULL) < 0)
        {
            cerr << "Error catching SIGALRM." << endl;
            exit(-1);
        }
    }

    if (pid > 0)
        m_setPIDCommands.erase(pid);

    if (m_setPIDCommands.size() == 0)
    {
        RestoreSignalHandlers();
    }

    if (pExitStatus != NULL)
        *pExitStatus = status;
    return pid;
}

bool CommandLoader::Wait(WaitInfo& waitInfo)
{
    waitInfo.m_pid = Wait(&waitInfo.m_status, waitInfo.m_pid,
                          waitInfo.m_timeOut);
    return(waitInfo.m_pid != -1 && waitInfo.m_status != -1 &&
           WIFEXITED(waitInfo.m_status) > 0);
}

CommandLoader* CommandLoader::GetInstance()
{
    if ( CommandLoader::s_pLoader == NULL )
        CommandLoader::s_pLoader = new CommandLoader();

    return CommandLoader::s_pLoader;
}

void CommandLoader::Destroy()
{
    delete CommandLoader::s_pLoader;
    s_pLoader = NULL;
}

void CommandLoader::StopAllRuns(int signo)
{
    SignalAllRuns(signo);

    // wait for the children to exit
    while (m_setPIDCommands.size() > 0)
    {
        Wait();
    }
}

void CommandLoader::SignalAllRuns(int signo) const
{
    // signal all of our children to exit
    for (set<pid_t>::iterator it = m_setPIDCommands.begin();
         it != m_setPIDCommands.end(); ++it)
    {
        kill(*it, signo);
    }
}

void CommandLoader::Log(const String& strMessage) const
{
    write(STDOUT_FILENO, strMessage.c_str(), strMessage.size());
    fsync(STDOUT_FILENO);
}

AssembleCommandLoader::AssembleCommandLoader( const String &strLogFile )
    : CommandLoader(), m_pidTail(-1)
{
    // create a pipe for reading and writing
    int fd[2];
    if (pipe(fd) < 0)
    {
        cerr << "Unable to create pipe" << endl;
        exit(0);
    }
    m_fdRead = fd[0];
    m_fdWrite = fd[1];

    // launch the tail process
    m_pidTail = LaunchTail(strLogFile);
}


pid_t AssembleCommandLoader::LaunchTail( const String &strLogFile )
{
    String strAssembleLoggerPath;
    if ( IsRegularFile( "./AssembleLogger" ) )
        strAssembleLoggerPath = "./AssembleLogger";
    else
        strAssembleLoggerPath = search_path_for_command("AssembleLogger");

    if ( strAssembleLoggerPath.empty() )
    {
        cout << "I cannot find AssembleLogger and I need it.  Exiting."
             << endl;
        exit(-1);
    }

    pid_t pid = fork();

    if (pid == -1)
    {
        cerr << "Unable to fork new process.  Exiting..." << endl;
        exit(0);
    }

    if (pid != 0)
    {
        // close the reading end of the pipe because the parent
        // has no use for it
        close(m_fdRead);
        return pid;
    }

    // close the write end of the pipe because the child have no use for it
    close(m_fdWrite);

    // swap out stdin
    dup2(m_fdRead, STDIN_FILENO);
    execl(strAssembleLoggerPath.c_str(), strAssembleLoggerPath.c_str(), strLogFile.c_str(),
          (char*) 0);
    exit(0);
}

bool AssembleCommandLoader::ShutDown()
{
    StopAllRuns();

    if (m_pidTail != -1)
    {
        //        kill(m_pidTail, SIGUSR1);
        close(m_fdWrite);

        // we'll wait until the tailing process exits to ensure it has
        // handled everything from the pipe
        while (waitpid(m_pidTail, NULL, 0) == -1 && errno == EINTR)
            continue;
        m_pidTail = -1;

    }

    return true;
}

void AssembleCommandLoader::Log(const String& strMessage) const
{
    write(m_fdWrite, strMessage.c_str(), strMessage.size());
}
