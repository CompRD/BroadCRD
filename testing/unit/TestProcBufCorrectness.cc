// A test of the speed of procbuf.

#include "MainTools.h"
#include "TaskTimer.h"

int main( int argc, char** argv )
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String( FILE );
  CommandArgument_UnsignedInt_OrDefault( ALLOC_MB, 0 );
  EndCommandArguments;


  size_t alloc_size = ALLOC_MB;
  alloc_size *= 1024 * 1024;
  PRINT( alloc_size );

  longlong *allocation = (longlong*) ( ALLOC_MB > 0 ? malloc( alloc_size ) : 0 );
  longlong *ptr = allocation;

  for ( size_t alloc_clear = 0; alloc_clear < alloc_size; alloc_clear += sizeof(longlong), ++ptr )
    *ptr = 0;

  const int bufferSize = 1024;
  char *fileBuffer = new char[ bufferSize ];
  char *procBuffer = new char[ bufferSize ];

  bool fail = false;
  
  ifstream fileIn( FILE.c_str() );

  String strCommand( "cat " + FILE );
  procbuf *pProcBuf = new procbuf( strCommand.c_str(), ios::in );
  istream procIn( pProcBuf );

  while ( fileIn && procIn ) 
  {
    fileIn.getline( fileBuffer, bufferSize );
    procIn.getline( procBuffer, bufferSize );
    if ( strcmp( fileBuffer, procBuffer ) != 0 )
    {
      printf( "difference found:\n" );
      printf( "from ifstream:        %s\n", fileBuffer );
      printf( "from procbuf/istream: %s\n", procBuffer );
      
      fail = true;
    }
  }

  if ( fileIn )
  {
    cout << "ifstream still good." << endl;
    fail = true;
  }
  if ( procIn )
  {
    cout << "procbuf/istream still good." << endl;
    fail = true;
  }

  if ( ! fail )
    cout << "Passed." << endl;
  else
    cout << "Failed." << endl;
  
  fileIn.close();
  delete pProcBuf;
}
