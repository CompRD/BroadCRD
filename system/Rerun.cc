/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// Rerun.

const char* DOC =

"Read output of command to find PrintTheCommandPretty header, reconstruct "
"the command, execute it."

;

#include "FastIfstream.h"
#include "MainTools.h"

void NoCommand(const char *msg)
{    
    cout << "Failed to find command.\n" << "Aborted with message: " << msg << "\n\n";
    exit(1);
}

int main(int argc, char *argv[]) 
{
    RunTime();

    BeginCommandArguments;
    CommandDoc(DOC);
    CommandArgument_String_OrDefault_Doc(F, "", "file to find command in");
    CommandArgument_String_OrDefault_Doc(TO, "", 
        "file to find command in and write output to");
    CommandArgument_String_OrDefault_Doc(C, "",
        "name of command; if not specified, use first encountered");
    CommandArgument_Int_OrDefault_Doc(N, 1, "instance of command");
    CommandArgument_Bool_OrDefault_Doc(GO, True,
        "if False, print command instead of executing it");
    CommandArgument_Bool_OrDefault_Doc(FORCE, False,
        "if True, override innate caution resulting from indeterminacy of "
        "command parsing prior to 2/9/08");
    EndCommandArguments;

    ForceAssertGe( N, 1 );

    ForceAssert(F == "" ^ TO == "");
    if (TO != "") 
        ForceAssert(GO);
    String infile = (F != "" ? F : TO);
    fast_ifstream in(infile);

    String line;
    String dashes(80, '-');

restart:
    
    // Search for a line of dashes demarcating the start of a command header.
    while(1)
    {    
        getline(in, line);
        if (in.fail()) 
            NoCommand("EOF while searching for prologue demarcation.");
        if (line == dashes)
            break;
    }
    
    // Look for the next line, containing the command timestamp and other info.
    getline(in, line);
    if (in.fail()) 
        NoCommand("EOF while scanning after dash.");

    if ( !( line.Contains("pid=") && line.Contains(" run ") ) )
                NoCommand("No keywords found."); 

    bool old = line.Contains(", using");
    String commandx;
    
    // Look for the command name.
    while(1)
    {
        getline(in, line);
        if (in.fail()) 
            NoCommand("EOF while searching for command.");
        // If we are looking for a specific command, then we should scan for it.
        if (C != "" && !line.Contains(" ", 0) && !line.Contains(C + " ", 0))
        {    
            while(1)
            {    
                getline(in, line);
                if (in.fail()) 
                    NoCommand("EOF while searching for specific command.");
                if (line.StartsWith(dashes)) 
                    break;
            }
            // We didn't encouter the command we were looking for,
            // keep looking.
            goto restart;
        }
        while(line.Contains(" ", 0)) 
            line = line.After(" ");
        if (!line.Contains(" \\", -1))
        {   
            commandx += line;
            break;    
        }
        line = line.RevBefore(" \\");
        while (line.Contains("  ", -1)) 
            line = line.RevBefore(" ");
        bool blank_end = line.Contains(" ", -1);
        if (old && !blank_end && !FORCE)
        {
            cout << "Unable to definitely determine original command.\n";
            cout << "Use FORCE=True to override.\n";
            cout << "Abort.\n\n";
            exit(1);    
        }
        commandx += line;
    }
    
    // Look for a line of dashes demarcating the end of the command header.
    getline(in, line);
    if (in.fail()) 
        NoCommand("EOF while searching for postscript demarcation.");
    if (!line.StartsWith(dashes)) 
        NoCommand("Could not find postscript demarcation.");
    if ( N > 1 )
    {    N--;
         goto restart;    }
    if (!GO) 
        cout << commandx << endl;
    else 
    {
        if (F != "")
            SystemSucceed(commandx);
        else 
            SystemSucceed(commandx + " > " + TO + " 2>&1");
    }
}
