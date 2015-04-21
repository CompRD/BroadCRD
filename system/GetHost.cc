// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

// Define GetDataDir( ), which returns the path prefix for the data directory.

#include "system/GetHost.h"
#include "system/System.h"
#include <cstdlib>

String GetDataDir()
{
    char* arachne_data_dir = getenv("ARACHNE_PRE");
    if ( !arachne_data_dir )
        FatalErr("In order to use this program, you must set the environment\n"
                "variable ARACHNE_PRE to the directory where you keep\n"
                "your assembly data files.  Please see the manual for details\n"
                "on setting this up, and then type\n\n"
                "     setenv ARACHNE_PRE x\n\n"
                "where \"x\" is the full path of the assembly data directory.");
    return arachne_data_dir;
}
