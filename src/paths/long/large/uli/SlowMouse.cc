///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/long/MakeKmerStuff.h"

int main( )
{
     String dir = "/wga/scr4/jaffe/GapToy/mouse.1456abcde/a.final";

     cout << Date( ) << ": loading" << endl;
     vecbasevector edges( dir + "/a.fastb" );

     cout << Date( ) << ": making lookup" << endl;
     const int K = 60;
     vec< kmer<K> > kmers;
     {    vec< triple<kmer<K>,int,int> > kmers_plus;
          MakeKmerLookup2( edges, kmers_plus );
          kmers.resize( kmers_plus.size( ) );
          #pragma omp parallel for
          for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
               kmers[i] = kmers_plus[i].first;    }
     ParallelUniqueSort(kmers);
     vec< vec<int64_t> > count(8);
     for ( int i = 0; i < 8; i++ )
          count[i].resize( kmers.size( ), 0 );

     vec<String> dx = { "52278.mouse1", "52280.mouse4", "52278.mouse5",
          "52278.mouse6", "52283.mousea", "52278.mouseb", "52348.mousec",
          "52348.moused" };
     cout << Date( ) << ": counting" << endl;
     #pragma omp parallel for
     for ( int idx = 0; idx < dx.isize( ); idx++ )
     {    String d = dx[idx];
          cout << Date( ) << ": loading " << idx+1 << " of " << dx.size( ) << endl;
          vecbasevector bases( 
               "/wga/scr4/jaffe/GapToy/" + d + "/data/frag_reads_orig.fastb" );
          kmer<K> x, xrc;
          cout << Date( ) << ": counting " << idx+1 << endl;
          for ( int64_t j = 0; j < (int64_t) bases.size( ); j++ )
          {    const basevector& b = bases[j];
               for ( int i = 0; i <= b.isize( ) - K; i++ )
               {    x.SetToSubOf( b, i );
                    xrc = x;
                    xrc.ReverseComplement( );
                    Bool fw = ( x < xrc );
                    if ( !fw ) x = xrc;
                    int64_t p = BinPosition( kmers, x );
                    if ( p >= 0 ) count[idx][p]++;    }    }    }

     cout << Date( ) << ": adding counts" << endl;
     #pragma omp parallel for
     for ( int64_t i = 0; i < kmers.jsize( ); i++ )
          count[1][i] += count[4][i];

     cout << Date( ) << ": creating kmers_plus again" << endl;
     vec< triple<kmer<K>,int,int> > kmers_plus;
     MakeKmerLookup2( edges, kmers_plus );

     cout << Date( ) << ": looking for responder kmers" << endl;
     vec<int> resp = {1,2,6};
     vec<int> nonresp = {0,3,5,7};
     #pragma omp parallel for
     for ( int64_t i = 0; i < kmers.jsize( ); i++ )
     {    
          // Require kmer seen at least 8 times in each responder and that the
          // total responder count is at least 20 times the total nonresponder
          // count.

          int64_t min_resp = 1000000000;
          int64_t sum_resp = 0, sum_nonresp = 0;
          for ( int j = 0; j < resp.isize( ); j++ )
          {    min_resp = Min( min_resp, count[ resp[j] ][i] );
               sum_resp += count[ resp[j] ][i];    }
          for ( int j = 0; j < nonresp.isize( ); j++ )
               sum_nonresp += count[ nonresp[j] ][i];
          if ( min_resp < 8 ) continue;
          if ( sum_resp < 20 * sum_nonresp ) continue;

          int64_t low = LowerBound1( kmers_plus, kmers[i] );
          int64_t high = UpperBound1( kmers_plus, kmers[i] );
          for ( int64_t j = low; j < high; j++ )
          {    int e = kmers_plus[j].second;
               PRINT2( kmers[i].ToString( ), e );    }    }

     Scram(0);    }
