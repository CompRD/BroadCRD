#include "MainTools.h"
#include "TaskTimer.h"

#define PREALLOC_MIN_SIZE 10
#define PREALLOC_MAX_SIZE 70
#define PREALLOC_SCALE 1024
#define NUM_CHUNKS 5000000
#define CHUNK_SIZE 120
#define NUM_ROUNDS 10


int main( )
{
  RunTime();

  typedef unsigned int  uint;

  TaskTimer timer;

  uint* prealloc_array;
  uint*** dynamic_ptr_array = new uint**[NUM_ROUNDS];
  for (uint round = 0; round < NUM_ROUNDS; round++ )
    dynamic_ptr_array[round] = new uint*[NUM_CHUNKS];

  for (uint var= PREALLOC_MIN_SIZE;
       var < PREALLOC_MAX_SIZE;
       var++)
    {
      cout << "Preallocating an array of size " << var << endl;
      timer.Start();
      prealloc_array = new uint[var*PREALLOC_SCALE];
      timer.Stop();
      cout << timer << endl;
      timer.Reset();

      for (uint round = 0; round < 2; round++)
	{
	  cout << "Allocating " << NUM_CHUNKS
	       << " chunks of size at least " << CHUNK_SIZE
	       << ": round " << round+1 << endl;

	  for (uint chunk = 0; chunk < NUM_CHUNKS; chunk++)
	    {
	      timer.Start();
	      dynamic_ptr_array[round][chunk] = new uint[CHUNK_SIZE];
	      timer.Stop();
	    }
	  cout << timer << endl;
	  timer.Reset();
	}

      cout << "Doing " << NUM_ROUNDS - 2 << " allocs and 2 deallocs per chunk: " << endl;
      TaskTimer newTimers[NUM_ROUNDS - 2], deleteTimer;
	  for (uint chunk = 0; chunk < NUM_CHUNKS; chunk++)
	    {
	      for ( uint round = 2; round < NUM_ROUNDS; ++round )
		{
	          newTimers[round-2].Start();
	          dynamic_ptr_array[round][chunk] = new uint[2*CHUNK_SIZE/round];
	          newTimers[round-2].Stop();
                }
	      for ( uint round = 0; round < 2; ++round )
		{
		  deleteTimer.Start();
		  delete dynamic_ptr_array[round][chunk];
		  deleteTimer.Stop();
		}

	      if ( chunk % 100000 == 0 )
		PrintMemUsage();
	    }
          for ( uint round = 2; round < NUM_ROUNDS; ++round )
	    cout << "News[" << round-2 << "]:\t" << newTimers[round-2] << endl;
	  cout << "Deletes:\t" << deleteTimer << endl;

      cout << "Freeing up all remaining allocated memory" << endl;
      timer.Start();
      delete[] prealloc_array;
      timer.Stop();
      for (uint round = 2; round < NUM_ROUNDS; round++)
	{
	  for (uint chunk = 0; chunk < NUM_CHUNKS; chunk++)
	    {
	      timer.Start();
	      delete[] dynamic_ptr_array[round][chunk];
	      timer.Stop();
	    }
	}
      cout << timer << endl;
      timer.Reset();

      cout << "------------------------------------------------------------" << endl;
    }

  cout << endl << "Timing test completed." << endl;
}

