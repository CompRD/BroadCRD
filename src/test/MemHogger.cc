///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MemHogger
//
// Test if we can allocate 90% of available memoery before it crashes or killed by
// an OOM killer.
//

#include "MainTools.h"

int main(int argc, char *argv[])
{
    RunTime( );
    // Check command line arguments
    // Thread control, etc.
    double clock = WallClockTime( );
    ostream & log = cout;


    const size_t mem_total = GetMaxMemory();
    const size_t mem_used  = MemUsageBytes();
    const size_t mem_avail = MemAvailable();

    cout << Date() << ": mem_total     = " << setw(18) << ToStringAddCommas(mem_total) << endl;
    cout << Date() << ": mem_used      = " << setw(18) << ToStringAddCommas(mem_used) << endl;
    cout << Date() << ": mem_avail     = " << setw(18) << ToStringAddCommas(mem_avail) << endl;

    // Try to allocate up to 90% of total memory using many 1GB vectors
    cout << Date() << ": Try to allocate up to 90% of total memory. " << endl;
    size_t alloc_max = mem_total * 9 / 10;

    { 
        size_t mem_used_new = mem_used;
        size_t alloc_step = mem_total / 20;
        vec< char* > mem_hogger;
        while ( mem_used_new < alloc_max ) {
            mem_used_new += alloc_step;
            cout << Date() << ":     Allocating " << ToStringAddCommas ( mem_used_new ) 
                << " bytes ... " << flush;
            char* test = new char[ alloc_step ];
            // access sparsely the memory
            for( size_t i = 0; i < alloc_step; i += 10000 ) {
                test[i] = '0';
            }
            cout << " OK." << endl;
            if ( test == 0 ) {
                cout << "Allocation failed " << endl;
                exit(1);
            }
            mem_hogger.push_back( test );
        }
        for ( size_t i = 0; i < mem_hogger.size(); ++i ) 
            delete[] mem_hogger[i];
    } // Local memory should be freed

    // Try to allocate up to 90% of total memory using a single vector
    cout << Date() << ": Try to allocate up to 90% of total memory using a single vector " << endl;
    { 
        size_t mem_used_new = mem_used;
        size_t alloc_step = mem_total / 10;
        vec<char> mem_hogger;
        mem_hogger.reserve( alloc_max );
        while ( mem_used_new < alloc_max ) {
            mem_used_new += alloc_step;
            cout << Date() << ":     Allocating " << ToStringAddCommas ( mem_used_new ) 
                << " bytes ... " << flush;
            mem_hogger.resize( mem_used_new );
            cout << " OK." << endl;
        }
    }

    cout << "\n" << Date( ) << ": time used = " << TimeSince(clock) << endl;
    cout << Date() << ": MemHogger done!" << endl;    
}
