/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file Qltout2SAMUtils.cc
 * \author tsharpe
 * \date Jan 13, 2009
 *
 * \brief Some utility classes to facilitate SAM creation.
 */
#include "lookup/Qltout2SAMUtils.h"
#include "solexa/SolexaDirectory.h"
#include "system/ErrNo.h"
#include "system/file/SymLink.h"
#include <sys/wait.h>

std::string qltout2sam::NULL_STR("NULL");

void qltout2sam::linkRefDict( File const& refDict, File const& refFasta )
{
    SymLink refFastaLink(refFasta);
    if ( refFastaLink.isValid() )
    {
        File dict(refFastaLink.target().changeExtension(SolexaDir::DICT_TYPE));
        if ( dict.exists() )
        {
            SymLink(refDict).setTarget(dict);
        }
    }
}

bool qltout2sam::nextAlignment( istream& is, Alignment& alignment, uint nReads, uint nRefs )
{
    bool result = false;
    
    std::string buf;
    while ( !result && is.good() )
    {
        getline(is, buf);
        if ( buf.size() > 6 && !buf.compare(0, 6, "QUERY\t") )
        {
            istringstream iss(buf.substr(6));
            iss >> alignment;

            string reason;
            if ( iss.fail() || iss.bad() || !iss.eof() )
            {
                cout << "Skipping alignment due to parse failure: " << buf << endl;
            }
            else if ( !alignment.isValid(reason,nReads,nRefs) )
            {
                cout << "Skipping alignment due to invalid field relationships (" << reason << "): " << buf << endl;
            }
            else
            {
                result = true;
            }
        }
    }
    return result;
}

using qltout2sam::Executor;

vec<Executor*> gChildren;

void qltout2sam::execChild( string command )
{
    gChildren.push_back(new Executor(command));
}

void qltout2sam::signalCatcher( int sigNo, siginfo_t* info, void* context )
{
    // ignore SIGCHLDs unless they indicate a "bad" child process termination
    if (sigNo == SIGCHLD &&
        ((info->si_code == CLD_EXITED && info->si_status == 0) ||
        (info->si_code != CLD_EXITED && info->si_code != CLD_KILLED 
        && info->si_code != CLD_DUMPED)))
    {
        return;
    }

    vec<Executor*>::iterator end = gChildren.end();
    for ( vec<Executor*>::iterator itr = gChildren.begin(); itr != end; ++itr )
    {
        Executor* pChild = *itr;
        if ( pChild )
        {
            if ( sigNo )
            {
                pChild->kill(sigNo);
            }
            delete pChild;
            *itr = 0;
        }
    }
    if ( sigNo != SIGUSR1 && sigNo != SIGUSR2 )
    {
        std::cerr << strsignal(sigNo) << std::endl;
        std::cerr << "Exiting" << std::endl;
        exit(sigNo);
    }
}

void qltout2sam::waitForChildren( string exitMessage )
{
    while ( true )
    {
        pid_t pid;
        int status;
        while ( (pid = wait(&status)) == -1 )
        {
            if ( errno == ECHILD ) // no more children to wait for
            {
                vec<Executor*>::iterator end = gChildren.end();
                for ( vec<Executor*>::iterator itr = gChildren.begin(); itr != end; ++itr )
                {
                    delete *itr;
                }
                return;
            }
            else if ( errno != EINTR )
            {
                ErrNo err;
                cout << "Wait call returned error: " << err << endl;
                cout << exitMessage << endl;
                signalCatcher(SIGTERM, NULL, NULL);
            }
        }

        vec<Executor*>::iterator end = gChildren.end();
        for ( vec<Executor*>::iterator itr = gChildren.begin(); itr != end; ++itr )
        {
            Executor* pChild = *itr;
            if ( pChild->reportStatus(pid,status) ) // if child terminated abnormally
            {
                cout << exitMessage << endl;
                signalCatcher(SIGTERM, NULL, NULL);
            }
        }
    }
}

uint qltout2sam::Read::gUnmappedRefID = INT_MAX;
int qltout2sam::NamedPipe::gNPipes;
char const* qltout2sam::Executor::EXEC_ARGS[] = { "/bin/sh", "-c", "command", 0 };

