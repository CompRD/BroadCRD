///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// JustSortKmers.  For 60x coverage of the human genome by 10 kb reads having 10% 
// substitution rate, find and sort all the 24-mers.  Report computational 
// performance.  Hardcoded to use 64 passes.
//
// In principle this could be used to determine the multiplicity of each kmer.
//
// If you switch to K=60 the total run time on crd9 is 4.53 hours.

#define _GLIBCXX_PARALLEL

#include <omp.h>

#include "Basevector.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "kmers/KmerRecord.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/DataSpec.h"
#include "paths/long/Logging.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/MakeSimReads.h"
#include "paths/long/Simulation.h"
#include "paths/long/local/Setup.h"

inline void GetPassIds( const basevector& u, const int L, const int j, 
     const int npasses, int& pass1, int& pass2 )
{    if ( npasses == 4 )
     {    pass1 = u[j];
          pass2 = 3 - u[j+L-1];    }
     else if ( npasses == 16 )
     {    pass1 = 4 * u[j] + u[j+1];
          pass2 = 4 * ( 3 - u[j+L-1] ) + ( 3 - u[j+L-2] );    }
     else if ( npasses == 64 )
     {    pass1 = 16 * u[j] + 4 * u[j+1] + u[j+2];
          pass2 = 16 * ( 3 - u[j+L-1] ) + 4 * ( 3 - u[j+L-2] )
               + ( 3 - u[j+L-3] );    }
     else if ( npasses == 256 )
     {    pass1 = 64 * u[j] + 16 * u[j+1] + 4 * u[j+2] + u[j+3];
          pass2 = 64 * ( 3 - u[j+L-1] ) + 16 * ( 3 - u[j+L-2] )
               + 4 * ( 3 - u[j+L-3] ) + ( 3 - u[j+L-4] );    }
     else
     {    FatalErr("npasses illegal");    }    }

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault_Doc(X, "", 
          "region definition, as in LongProto");
     EndCommandArguments;

     // Initial set up.

     // Create simulated reads.  Rather ugly.

     long_sim sim( "{LEN=10000,ERR_DEL=0,ERR_INS=0,ERR_SUB=0.1}" );
     vecbasevector bases;
     {    long_data_spec spec( "" );
          String X_actual;
          long_logging logc( "", "" );
          RefTraceControl RTCtrl( "" );
          ref_data ref;
          BuildRefData( "", "human", X, spec.HUMAN_CONTROLS, X_actual, ref,
               0, logc, RTCtrl );
          vec<ref_loc> readlocs;
          MakeSimReads( sim, ref, bases, readlocs );    }

     // Go through passes.

     double clock = WallClockTime( );
     const int npasses = 64;
     for ( int pass = 0; pass < npasses; pass++ )
     {    cout << "starting pass " << pass+1 << " of " << npasses << endl;

          // Index kmers in the reads.

          const int K = 24;
          vec< kmer<K> > kmers;
          vec<int> counts( bases.size( ), 0 );
          #pragma omp parallel for
          for ( size_t i = 0; i < bases.size( ); i++ )
          {    const basevector& u = bases[i];
               int count = 0;
               for ( int j = 0; j <= u.isize( ) - K; j++ )
               {    int pass1 = 0, pass2 = 0;
                    GetPassIds( u, K, j, npasses, pass1, pass2 );
                    if ( pass1 < pass || pass2 < pass ) continue;
                    if ( pass1 > pass && pass2 > pass ) continue;
                    counts[i]++;    }    }
          vec<int64_t> starts( bases.size( ) + 1 );
          starts[0] = 0;
          for ( size_t i = 0; i < bases.size( ); i++ )
               starts[i+1] = starts[i] + counts[i];
          kmers.resize( starts.back( ) );
          #pragma omp parallel for
          for ( size_t i = 0; i < bases.size( ); i++ )
          {    const basevector& u = bases[i];
               int count = 0;
               kmer<K> x, xrc;
               for ( int j = 0; j <= u.isize( ) - K; j++ )
               {    int pass1 = 0, pass2 = 0;
                    GetPassIds( u, K, j, npasses, pass1, pass2 );
                    if ( pass1 < pass || pass2 < pass ) continue;
                    if ( pass1 > pass && pass2 > pass ) continue;
                    Bool type1 = ( pass1 == pass ), type2 = ( pass2 == pass );
                    int64_t r = starts[i] + count++;
                    x.SetToSubOf( u, j );
                    xrc = x;
                    xrc.ReverseComplement( );   
                    if ( x <= xrc ) kmers[r] = x;
                    else kmers[r] = xrc;    }    }
          ParallelSort(kmers);    }
     cout << Date( ) << ": done, time used = " << TimeSince(clock)
          << ", peak mem used = " << setiosflags(ios::fixed) << setprecision(1)
          << PeakMemUsageBytes( ) / 1000000000.0 << resetiosflags(ios::fixed)
          << " GB" << endl;    }
