/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// EffectOfGcWindowBias.  Suppose given a genome and the distribution of the GC
// contents of the 50-base windows beginning at read start points. Suppose reads are
// generated at random relative to this constraint.  Determine how biased the 
// resulting data set would be.
//
// GENOME = fastb file for genome
// GC = {n0,...,n50} were ni = percent of start positions for which number of GCs
// in the 50-base window is i.
// SAMPLE = sample size = number of reads generated.

#include "Basevector.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "bias/UniformBias.h"
#include "math/Functions.h"
#include "random/Random.h"

int main( int argc, char *argv[] ) 
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(GENOME);
     CommandArgument_String(GC);
     CommandArgument_Int(SAMPLE);
     EndCommandArguments;

     vec<double> gc;
     ParseDoubleSet( GC, gc, False );
     ForceAssertEq( gc.size( ), 51u );
     for ( int i = 0; i < 51; i++ )
          gc[i] /= 100.0;

     vecbasevector genome(GENOME);
     vecbasevector genome_rc(genome);
     ReverseComplement(genome_rc);
     genome.Append(genome_rc);

     vec<unsigned int> sum( genome.size( ) + 1 );
     sum[0] = 0;
     for ( size_t i = 0; i < genome.size( ); i++ )
          sum[i+1] = sum[i] + genome[i].size( );
     int N = sum[ genome.size( ) ];

     vec<int> counts( 51, 0 );
     vec< vec<int> > starts( genome.size( ) );
     for ( size_t i = 0; i < genome.size( ); i++ )
          starts[i].resize( genome[i].size( ), 0 );
     for ( int i = 0; i < SAMPLE; i++ )
     {    int d = 0;
          if ( i == 0 ) d = Position( gc, Max(gc) );
          else
          {    static vec<double> cd(51);
               for ( int j = 0; j < 51; j++ )
                    cd[j] = double( counts[j] ) / double(i);
               for ( int j = 1; j < 51; j++ )
                    if ( gc[j] - cd[j] > gc[d] - cd[d] ) d = j;    }
          unsigned int tig, n;
          while(1)
          {    n = randomx( ) % N;
               for ( tig = 0; tig < genome.size( ); tig++ )
                    if ( sum[tig] > n ) break;
               tig--;
               n -= sum[tig];
               if ( n+50 > genome[tig].size( ) ) continue;
               int x = 0;
               for ( int k = 0; k < 50; k++ )
                    if ( genome[tig][n+k] == 1 || genome[tig][n+k] == 2 ) ++x;
               if ( x == d ) break;    }
          ++starts[tig][n];
          ++counts[d];    }

     vec<int> sizes;
     for ( size_t i = 0; i < genome.size( ); i++ )
          sizes.push_back( genome[i].size( ) );
     cout << HowBiased( starts, sizes ) << " +/- "
          << HowBiasedDev( starts, sizes ) << "\n";    }
