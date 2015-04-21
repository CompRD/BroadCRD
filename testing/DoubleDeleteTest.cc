#include <signal.h>

#include "system/MemTracker.h"

#include "MainTools.h"


void simple_arachne_signal_handler( int signal_number, siginfo_t* info,
    void* context )
{
  if ( signal_number == SIGTERM )
      _exit(1);

  cout << "\nNot generating a backtrace..." << endl;
  _exit(1);
}


int main(int argc, char** argv )
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String( handler );
  CommandArgument_Bool_OrDefault( double_delete, True );
  CommandArgument_Bool_OrDefault( missing_delete, True );
  CommandArgument_Bool_OrDefault( array_delete_object, True );
  CommandArgument_Bool_OrDefault( object_delete_array, True );
  EndCommandArguments;

  ArachneSignalHandler *pHandlerFunc = 0;
  if ( handler == "normal" )
    pHandlerFunc = &arachne_signal_handler_standard;
  else if ( handler == "simple" )
    pHandlerFunc = &simple_arachne_signal_handler;
  else if ( handler == "none" )
    pHandlerFunc = 0;
  else
  {
    cout << "The handler option should be one of 'normal', 'simple', or 'none'." << endl;
    exit(0);
  }

  ArachneInterruptHandler( pHandlerFunc );
  
  const int numPtrs = 10;
  int **x = new int* [numPtrs];
  for ( int ptrNum = 0; ptrNum < numPtrs; ++ptrNum )
    x[ptrNum] = new int;

  for ( int ptrNum = 1; ptrNum < numPtrs; ++ptrNum )
    delete x[ptrNum];

  if ( double_delete )
    delete x[1];

  if ( ! missing_delete )
    delete x[0];

  if ( object_delete_array )
  {
    int *y = new int[10];
    delete y;
  }

  if ( array_delete_object )
  {
    int *y = new int;
    delete [] y;
  }

  delete [] x;

  return 0;
}

