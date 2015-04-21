// A test of the speed of procbuf.

#include "MainTools.h"
#include "TaskTimer.h"

int main( int argc, char** argv )
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String( FILE );
  EndCommandArguments;

  
  const int bufferSize = 1024;
  char *buffer = new char[ bufferSize ];
  
  TaskTimer timer;

  ifstream fileIn( FILE.c_str() );

  PRINT( Date() );

  timer.Start();
  while ( fileIn ) 
    fileIn.getline( buffer, bufferSize );
  timer.Stop();

  fileIn.close();
  
  cout << "ifstream:\t\t" << timer << endl;
  timer.Reset();

  String strCommand( "cat " + FILE );
  procbuf *pProcBuf = new procbuf( strCommand.c_str(), ios::in );
  istream procIn( pProcBuf );
  
  timer.Start();
  while ( procIn )
    procIn.getline( buffer, bufferSize );
  timer.Stop();
  
  cout << "procbuf/istream:\t" << timer << endl;

  delete pProcBuf;
}
