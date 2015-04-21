// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

// DupCount.  Draw n items from a set of N, replacing each time.  Compute the number
// of duplicates found (amongst the n).  Repeat this experiment "reps" times, 
// yielding a distribution.  Determine the placement of "test" in this distribution.

#include "MainTools.h"
#include "random/Random.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_UnsignedInt(N);
     CommandArgument_UnsignedInt(n);
     CommandArgument_UnsignedInt(reps);
     CommandArgument_UnsignedInt(test);
     EndCommandArguments;

     vec<unsigned int> dup_counts;
     for ( unsigned int p = 0; p < reps; p++ )
     {    vec<int> count( N, 0 );
          for ( unsigned int i = 0; i < n; i++ )
               ++count[ randomx( ) % N ];
          int dups = 0;
          for ( unsigned int i = 0; i < N; i++ )
               if ( count[i] > 1 ) dups += count[i];
          dup_counts.push_back(dups);    }
     Sort(dup_counts);
     for ( unsigned int i = 0; i < reps; i++ )
     {    if ( test <= dup_counts[i] )
          {    cout << PERCENT_RATIO( 3, i, reps ) << "\n";
               break;    }    }    }
