// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
//


#include "String.h"
#include "system/RunTime.h"
#include <iostream>

int main( int argc, char** argv )
{
  RunTime();

  String command_line, command;
  while ( 1 )
  {
    getline( cin, command_line );

    if ( ! cin )
      break;

    // append this line to command
    command += command_line; 

    if ( command_line.Contains( "\\", -1 ) )
    {
      // chop off backslash
      command.resize( command.size() - 1 ); 
      // get next line
      continue;
    }
  
    cout << system( command.c_str() ) << endl;

    command.resize(0);
  }

  exit(0);
}
