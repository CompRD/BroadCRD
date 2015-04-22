///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/**
   Program: FileIOTest
   
   Designed to determine the read and write throughput of a filesystem.
   It writes a file of the requested size and then tries read it back in.
   The file contents are streamed in and out, nothing is kept in memory.

   With the THREADS option and using OpenMP it is possible to split the test
   over multiple threads and files. Each thread will write to or read from
   a seperate file - the sum of the files equaling the requested test size.
   
   The threaded test allows you to determine if there are any limitations or
   bottlenecks around the number open files being accessed concurrently.
*/

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
#define _GLIBCXX_PARALLEL

#include <omp.h>

#include "MainTools.h"
#include "TaskTimer.h"
#include "system/HostName.h"

int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandDoc(
    "File system IO test with optional multithreaded IO test");
  CommandArgument_Bool_OrDefault_Doc(WRITE, True,
    "Perform write IO test");
  CommandArgument_Bool_OrDefault_Doc(READ, True,
    "Perform read IO test");
  CommandArgument_String_OrDefault_Doc(BASENAME, "speedtest",
    "Basename for test files");
  CommandArgument_Int_OrDefault_Valid_Doc(SIZE, 10, "[1,)",
    "Total IO test size in GB");
  CommandArgument_Int_OrDefault_Valid_Doc(THREADS, 1, "[1,)",
    "Number of threads (and files) over which to spread the test");
  EndCommandArguments;

  THREADS = configNumThreads(THREADS);
  omp_set_num_threads(THREADS);

  cout << Date() << " Starting I/O Tests" << endl << endl;

  // Compute number of chunks per required per thread
  const longlong chunk_size = 1024;
  const longlong file_size = SIZE;
  const longlong chunks = file_size * 1024 * 1024 * 1024 / (chunk_size * THREADS);


  // Write Test
  if (WRITE) {

    cout << Date() << " Writing " << SIZE << " GB of data using " << THREADS << " threads" << endl;
    TaskTimer timer;
    timer.Start();

    // thread local variables
    ofstream file_out;
    String filename;
    int tid;
    char* buffer;

    cout << "  Starting threads: ";

    #pragma omp parallel private(tid, file_out, filename, buffer)
    {
      tid = omp_get_thread_num();
      filename = BASENAME + +"_" + getHostName() + "_" + ToString(SIZE) + "GB_T" + ToString(tid);

      #pragma omp critical
      cout << tid << " " << flush;

      // Prepare individual buffer for file I/O operations
      buffer = new char[chunk_size];

      // Write to files
      file_out.open( filename.c_str( ), ios::out);
      for (long long i = 0; i < chunks; ++i) 
 	file_out.write(buffer, chunk_size);
      file_out.close();

      delete[] buffer;
    }

    timer.Stop();
    const float runtime = timer.Elapsed();

    cout << endl << Date() << " Write complete " << endl;
    cout << "  Runtime   : " << runtime << " secs" << endl;
    cout << "  I/O speed : " << SIZE * 1000 / runtime << " MB/sec" << endl << endl;

  }


  // Read Test
  if (READ) {
    cout << Date() << " Reading " << SIZE << " GB of data using " << THREADS << " threads" << endl;
    TaskTimer timer;
    timer.Start();

    // thread local variables
    ifstream file_in;
    String filename;
    int tid;
    char* buffer;

    cout << "  Starting threads: ";

    #pragma omp parallel private(tid, file_in, filename, buffer)
    {
      tid = omp_get_thread_num();
      filename = BASENAME + +"_" + getHostName() + "_" + ToString(SIZE) + "GB_T" + ToString(tid);

      #pragma omp critical
      cout << tid << " " << flush;

      // Prepare individual buffer for file I/O operations
      buffer = new char[chunk_size];

      // Read from files
      file_in.open( filename.c_str( ), ios::in);
      for (long long i = 0; i < chunks; ++i) 
	file_in.read(buffer, chunk_size);
      file_in.close();

      delete[] buffer;
    }

    timer.Stop();
    const float runtime = timer.Elapsed();

    cout << endl << Date() << " Read complete " << endl;
    cout << "  Runtime   : " << runtime << " secs" << endl;
    cout << "  I/O speed : " << SIZE * 1000 / runtime << " MB/sec" << endl << endl;
  }


  cout << Date() << " Testing Complete " << endl;
   
 }
