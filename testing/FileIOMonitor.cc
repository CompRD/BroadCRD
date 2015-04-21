///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/**
   Program: FileIOMonitor
   
   Module to monitor the read and write throughput of a filesystem.
   It just writes (or reads) a file of fixed size repeatedly.
*/



#include "MainTools.h"
#include "TaskTimer.h"

int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandDoc(
    "Very basic file system IO monitor");
  CommandArgument_Bool_OrDefault_Doc(READ, True,
    "Perform read IO test (otherwise perform write IO test)");
  CommandArgument_String_OrDefault_Doc(FILENAME, "monitor",
    "Basename for test files");
  CommandArgument_Int_OrDefault_Valid_Doc(SIZE, 500, "[1,)",
    "Total IO test size in MB");
  EndCommandArguments;

  cout << Date() << " Starting I/O Monitor" << endl;

  // Compute number of chunks per file
  const int chunk_size = 1024;  // Bytes
  const long long chunks = SIZE * 1024 * 1024 / chunk_size;

  TaskTimer timer;

  // Prepare buffer for file I/O operations
  char* buffer = new char[chunk_size];

  // Write Test (or prepare file for read test) 
  while (true) {
      
    ofstream file_out;
    file_out.open( FILENAME.c_str( ), ios::out);
    
    timer.Start();
    for (long long i = 0; i < chunks; ++i) 
      file_out.write(buffer, chunk_size);
    file_out.close();
    
    timer.Stop();
    const float runtime = timer.Elapsed();
    cout << Date() << ": " << SIZE * 60.0 / (runtime * 1024)<< " GB/min (write)" << endl;    
    timer.Reset();

    // For read test, only write file once and continue
    if (READ)
      break;
   
  }
  

  // Read Test
  if (READ) {

    while (true) {
      
      ifstream file_in;
      file_in.open( FILENAME.c_str( ), ios::in);

      timer.Start();
      for (long long i = 0; i < chunks; ++i) 
	file_in.read(buffer, chunk_size);
      file_in.close();

      timer.Stop();
      const float runtime = timer.Elapsed();
      cout << Date() << ": " << SIZE * 60.0 / (runtime * 1024)<< " GB/min (read)" << endl;
      timer.Reset();
	 
    }
  }

  cout << Date() << " Testing Complete " << endl;

 }
