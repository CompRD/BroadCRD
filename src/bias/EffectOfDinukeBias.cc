/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// EffectOfDinukeBias.  Suppose given a genome and the distribution of the
// dinucleotide pairs at the juncture of read start points (one base before
// the start and one after).  Suppose reads are generated at random relative
// to this constraint.  Determine how biased the resulting data set would be.
//
// GENOME = fastb file for genome
// DINUKE = {AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT} where each
//          nucleotide pair is replaced by a floating-point number representing
//          the percent frequency of that pair
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
     CommandArgument_String(DINUKE);
     CommandArgument_Int(SAMPLE);
     EndCommandArguments;

     vec<double> dinuke;
     ParseDoubleSet( DINUKE, dinuke, False );
     ForceAssertEq( dinuke.size( ), 16u );
     for ( int i = 0; i < 16; i++ )
          dinuke[i] /= 100.0;

     vecbasevector genome(GENOME);
     vecbasevector genome_rc(genome);
     ReverseComplement(genome_rc);
     genome.Append(genome_rc);

     vec<unsigned int> sum( genome.size( ) + 1 );
     sum[0] = 0;
     for ( size_t i = 0; i < genome.size( ); i++ )
          sum[i+1] = sum[i] + genome[i].size( );
     unsigned int N = sum[ genome.size( ) ];

     vec<int> counts( 16, 0 );
     vec< vec<int> > starts( genome.size( ) );
     for ( size_t i = 0; i < genome.size( ); i++ )
          starts[i].resize( genome[i].size( ), 0 );
     for ( int i = 0; i < SAMPLE; i++ )
     {    int d = 0;
          if ( i == 0 ) d = Position( dinuke, Max(dinuke) );
          else
          {    static vec<double> cd(16);
               for ( int j = 0; j < 16; j++ )
                    cd[j] = double( counts[j] ) / double(i);
               for ( int j = 1; j < 16; j++ )
                    if ( dinuke[j] - cd[j] > dinuke[d] - cd[d] ) d = j;    }
          unsigned int tig, n;
          while(1)
          {    n = randomx( ) % N;
               for ( tig = 0; tig < genome.size( ); tig++ )
                    if ( sum[tig] > n ) break;
               tig--;
               n -= sum[tig];
               if ( n == 0 ) continue;
               int base_before = genome[tig][n-1], base_after = genome[tig][n];
               if ( d == 4 * base_before + base_after ) break;    }
          ++starts[tig][n];
          ++counts[d];    }

     vec<int> sizes;
     for ( size_t i = 0; i < genome.size( ); i++ )
          sizes.push_back( genome[i].size( ) );
     cout << HowBiased( starts, sizes ) << " +/- "
          << HowBiasedDev( starts, sizes ) << "\n";    }
