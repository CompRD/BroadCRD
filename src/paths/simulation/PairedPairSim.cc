// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

// PairedPairSim: create random read pairs from specified libraries, attempt to
// walk across inserts using "paired pair" rules, compute the fraction of inserts
// which close.
//
// The paired pair rules work the following way.  You start with a read pair.
// Then you gradually add read pairs.  A read pair can be added if both ends
// overlap existing read pair bases by at least k.
//
// Arguments:
// G = genome size
// n = read length
// N = mean insert length
// d = standard deviation of insert length
// k = minimum overlap
// C = coverage
// s = number of seeds to try

#include "math/Functions.h"
#include "math/HoInterval.h"
#include "MainTools.h"
#include "random/NormalRandom.h"
#include "random/Random.h"

class read_datum {

     public:

     int start;
     int insert_id;
     Bool rc;

     read_datum( ) { }
     read_datum( int start, int insert_id, Bool rc )
          : start(start), insert_id(insert_id), rc(rc) { }

     friend Bool operator<( const read_datum& r1, const read_datum& r2 )
     {    return r1.start < r2.start;    }

};

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_Int(G);
     CommandArgument_Double(C);
     CommandArgument_Int(n);
     CommandArgument_Int(N);
     CommandArgument_Int(K);
     CommandArgument_Double(d);
     CommandArgument_Int(s);
     EndCommandArguments;

     // Create inserts.

     int npairs = int( floor( (double) G * C ) ) / (2*n);
     vec<ho_interval> inserts(npairs);
     vec<read_datum> reads( npairs*2 );
     NormalRandom ilength( N, d );
     for ( int i = 0; i < npairs; i++ )
     {    int len = int( round( ilength.value( ) ) );
          ForceAssertGt( len, 0 );
          ForceAssertLt( len, G );
          int start = randomx( ) % ( G - len );
          inserts[i] = ho_interval( start, start + len );    
          reads[ 2*i ] = read_datum( start, i, False );
          reads[ 2*i + 1 ] = read_datum( start + len - n, i, True );    }
     Sort(reads);

     // Build paired pair constructs.

     vec<Bool> tried( inserts.size( ), False );
     int closed = 0;
     for ( int i = 0; i < s; i++ )
     {    int pi = randomx( ) % npairs;
          if ( tried[pi] )
          {    --i;
               continue;    }
          tried[pi] = True;
          const ho_interval& I = inserts[pi];

          // Don't use inserts too near to the ends of the genome.

          int avoid = 10 * n;
          if ( I.Start( ) < avoid || G - I.Stop( ) < avoid ) 
          {    --i;
               continue;    }

          // There are two clusters of reads, represented as intervals, which can
          // grow.  Initially they are the end reads of the seed insert I.

          ho_interval left( I.Start( ), I.Start( ) + n );
          ho_interval right( I.Stop( ) - n, I.Stop( ) );

          // Grow the clusters.

          while(1)
          {    Bool progress = False;

               // Check to see if we're done.

               if ( left.Stop( ) - right.Start( ) >= K ) break;

               // Go through all the pairs which could extend the cluster.

               int low = lower_bound( reads.begin( ), reads.end( ),
                    read_datum( left.Start( ) - ( n - K ), 0, False ) )
                    - reads.begin( );
               int high = upper_bound( reads.begin( ), reads.end( ),
                    read_datum( right.Stop( ) - K, 0, False ) ) - reads.begin( );
               for ( int j = low; j < high; j++ )
               {    const ho_interval& m = inserts[ reads[j].insert_id ];
                    ho_interval mleft( m.Start( ), m.Start( ) + n );
                    ho_interval mright( m.Stop( ) - n, m.Stop( ) );
                    Bool own_mleft = Subset(mleft, left) || Subset(mleft, right);
                    Bool own_mright = Subset(mright, left) || Subset(mright, right);
                    if ( own_mleft && own_mright ) continue;
                    Bool mleft_adds = Overlap(mleft, left) >= K
                         || Overlap(mleft, right) >= K;
                    Bool mright_adds = Overlap(mright, left) >= K
                         || Overlap(mright, right) >= K;
                    if ( mleft_adds && mright_adds )
                    {    if ( Overlap(mleft, left) >= K ) 
                              left = Span(mleft, left);
                         if ( Overlap(mleft, right) >= K ) 
                              right = Span(mleft, right);
                         if ( Overlap(mright, left) >= K ) 
                              left = Span(mright, left);
                         if ( Overlap(mright, right) >= K ) 
                              right = Span(mright, right);    }

                    else if ( Overlap(mright, left) >= K 
                         && mright.Start( ) < left.Start( ) )
                    {    left.Set( mright.Start( ), left.Stop( ) );    }
                    else if ( Overlap(mleft, right) >= K
                         && mleft.Stop( ) > right.Stop( ) )
                    {    right.Set( right.Start( ), mleft.Stop( ) );    }

                    else continue;
                    progress = True;    }
               if ( !progress ) break;    }

          // Record results.

          if ( left.Stop( ) - right.Start( ) >= K ) ++closed;    
          else cout << "see gap of size " << right.Start( ) - left.Stop( ) << endl;

               }

     // Summarize.

     cout << PERCENT_RATIO( 4, closed, s ) << " closed\n";    }
