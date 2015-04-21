/** \file TestPthread.cc
An example program demonstrating the use of the Pthread class.
*/

#include "MainTools.h"
#include "system/Pthread.h"

/// A function with signature 'void *FunctionName(void *argument)' that will be called by the thread.  Pass a pointer to a struct or class to specify multiple arguments.
void *ThreadedFunction(void *threadid) {
    int *tid = (int *) threadid;

    int sleeptime = (rand() % 5) + 1;
    for (unsigned int i = 0; i < 3; i++) {
        sleep(sleeptime);
        cout << "This is thread number " << (*tid) << ", iteration " << i << "." << endl;
    }

    int *retval = new int;
    *retval = sleeptime;

    return retval;
}

int main(int argc, char **argv) {
    RunTime();

    BeginCommandArguments;
    CommandArgument_UnsignedInt(JOBS);
    EndCommandArguments;

    srand(time(NULL));

    // Instantiate the thread pool
    vecThread threads(JOBS);

    vec<unsigned int> ids(JOBS);
    for (unsigned int threadid = 0; threadid < threads.size(); threadid++) {
        ids[threadid] = threadid;

        // Bind a function and an argument to each thread
        threads[threadid].Bind(ThreadedFunction, &ids[threadid]);
    }

    // Execute all the threads
    threads.Run(JOBS/3);

    // Print out the return values from each thread.
    for (unsigned int threadid = 0; threadid < threads.size(); threadid++) {
        cout << "Thread " << threadid << " slept for " << *(reinterpret_cast<int *>(threads[threadid].GetReturnValue())) << " seconds." << endl;
    }

    return 0;
}
