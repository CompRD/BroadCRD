///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Simo.  Estimate error rates for consensus of simulated reads.
//
// See Simo.png for results.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <random>

#include <sys/time.h>

#include "MainTools.h"
#include "PrintAlignment.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/long/magic/BasicScore.h"
#include "paths/long/magic/IterCon.h"
#include "paths/long/ultra/ThreadedBlocks.h"
#include "random/Random.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;

     CommandArgument_Double_OrDefault(SUB, 21); // substitution error rate
     CommandArgument_Double_OrDefault(DEL, 12); // deletion error rate
     CommandArgument_Double_OrDefault(INS,  2); // insertion error rate
     CommandArgument_Double_OrDefault(COV, 12); // mean coverage

     CommandArgument_Double_OrDefault_Doc(MUL,  2,
          "used to define position-dependence of errors; smaller is worse");

     CommandArgument_Int_OrDefault(SAMPLE, 100);

     EndCommandArguments;

     vecbasevector genome( "/wga/dev/references/Escherichia_coli/genome.fastb" );

     int len = 100;

     vec<int> all_errs;
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int d = 0; d < SAMPLE; d++ )
     {
          timeval t;
          gettimeofday( &t, NULL );
          uint seed = t.tv_sec + t.tv_usec + d*1000000;

          vecbasevector segs;

          int start = randomx( ) % ( genome[0].isize( ) - len );

          basevector b( genome[0], start, len );
          basevector borig(b);

          std::default_random_engine generator;
          generator.seed(seed);
          std::poisson_distribution<int> distribution(COV);
     
          int covx = distribution(generator);

          poisson_distribution<int> muld(MUL);
          vec<double> mulx(len);
          for ( int i = 0; i < len; i++ )
               mulx[i] = muld(generator) / double(MUL);

          for ( int c = 0; c < covx; c++ )
          {    basevector x(b);
               for ( int i = b.isize( ) - 1; i >= 0; i-- )
               {    if ( randomx( ) % 10000 < 100*SUB*mulx[i] )
                         x.Set( i, ( x[i] + 1 + ( randomx( ) % 3 ) ) % 4 );
                    if ( randomx( ) % 10000 < 100*INS*mulx[i] )
                    {    basevector y;
                         for ( int j = 0; j < x.isize( ); j++ )
                         {    y.push_back( x[j] );
                              if ( j == i ) y.push_back( randomx( ) % 4 );    }
                         x = y;     }
                    if ( randomx( ) % 10000 < 100*DEL*mulx[i] )
                    {    basevector y;
                         for ( int j = 0; j < x.isize( ); j++ )
                              if ( j != i ) y.push_back( x[j] );
                         x = y;    }    }

               /*
               int best_loc;
               alignment a;
               if ( x.size( ) <= b.size( ) )
                    SmithWatFree( x, b, best_loc, a );
               else
               {    SmithWatFree( b, x, best_loc, a );
                    a.Flip( );    }
               PrintVisualAlignment( False, cout, x, b, a );
               */
               
               segs.push_back(x);    }

          basevector bx = IterCon(segs);
     
          int best_loc;
          alignment a;
          if ( bx.size( ) <= borig.size( ) )
               SmithWatFree( bx, borig, best_loc, a );
          else
          {    SmithWatFree( borig, bx, best_loc, a );
               a.Flip( );    }
          align al(a);
          #pragma omp critical
          {    // PRINT(covx);
               // cout << "\nnew alignment:\n";
               // PrintVisualAlignment( False, cout, bx, borig, a );
               // cout << "\n";
               int errs = al.Errors( bx, borig );
               // cout << "errs = " << errs << endl;    
               all_errs.push_back(errs);    }    }
     
     Sort(all_errs);
     vec<int> all_errs_sub;
     for ( int i = SAMPLE/4; i < 3*SAMPLE/4; i++ )
          all_errs_sub.push_back( all_errs[i] );
     cout << "mean errs = " << Mean(all_errs_sub) << endl;    }
