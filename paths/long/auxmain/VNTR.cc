///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// VNTR.  Given a MUC1-type region and given a seed repeat (hardcoded), attempt to 
// translate read pairs into concatenations of units, each having the size of
// the repeat.  Tested on one repeat suggested by Giulio Genovese.  Best results
// were obtained from dataset 3, but results were lousy, suggesting that the
// repeat interferes with cluster generation or sequencing, even though its
// sequence composition appears normal.

#include "Basevector.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "Qualvector.h"
#include "VecUtilities.h"

// LongProto SAMPLE=human READS=#picard TMP=tmp.xxx
// X=12:2364957-2365257 DATASET=3 EXIT=CORRECT

String Pic( const int x, const int n )
{    int d = ToString(n).size( );
     String p = ToString(x);
     return String( d - p.isize( ), '0' ) + p;    }

int main( )
{    RunTime( );

     vecbasevector bases( "tmp.xxx/frag_reads_mod0.fastb" );
     vecqualvector quals( "tmp.xxx/frag_reads_mod0.qualb" );
     PairsManager pairs;
     pairs.Read( "tmp.xxx/frag_reads_orig.pairs" );

     int N = bases.size( );

     basevector r( "CACACGACCCTGACCTGACTAGTTTACAAC" );
     int n = r.size( );

     vec<String> units;
     units.push_back( r.ToString( ) );
     int pos = 0;
     const int minc = 10;
     while(1)
     {    vec<vec<String>> left_neighbors( units.size( ) );
          vec<vec<String>> right_neighbors( units.size( ) );
          for ( int i = 0; i < N; i++ )
          for ( int pass = 1; pass <= 2; pass++ )
          {    basevector b = bases[i];
               if ( pass == 2 ) b.ReverseComplement( );
               String bs = b.ToString( );
               for ( int j = 0; j <= b.isize( ) - n; j++ )
               {    for ( int k = pos; k < units.isize( ); k++ )
                    {    if ( bs.Contains( units[k], j ) )
                         {    if ( j >= n )
                              {    left_neighbors[k].push_back( 
                                        bs.substr( j - n, n ) );    }
                              if ( j <= b.isize( ) - 2*n )
                              {    right_neighbors[k].push_back( bs.substr( 
                                        j + n, n ) );    }    }    }    }    }
          vec<String> neighbors;
          for ( int k = pos; k < units.isize( ); k++ )
          {    Sort( left_neighbors[k] );
               for ( int i = 0; i < left_neighbors[k].isize( ); i++ )
               {    int j = left_neighbors[k].NextDiff(i);
                    if ( j - i >= minc )
                    {    if ( !Member( units, left_neighbors[k][i] ) )
                              neighbors.push_back( left_neighbors[k][i] );    }
                    i = j - 1;    }
               Sort( right_neighbors[k] );
               for ( int i = 0; i < right_neighbors[k].isize( ); i++ )
               {    int j = right_neighbors[k].NextDiff(i);
                    if ( j - i >= minc )
                    {    if ( !Member( units, right_neighbors[k][i] ) )
                              neighbors.push_back( right_neighbors[k][i] );    }
                    i = j - 1;    }    }
          UniqueSort(neighbors);
          if ( neighbors.empty( ) ) break;
          pos = units.size( );
          units.append(neighbors);    }

     const int minq = 30;
     Sort(units);
     vec<int> counts( units.size( ), 0 );
     for ( int i = 0; i < N; i++ )
     for ( int pass = 1; pass <= 2; pass++ )
     {    basevector b = bases[i];
          qualvector q = quals[i];
          if ( pass == 2 ) 
          {    b.ReverseComplement( );
               q.ReverseMe( );    }
          String bs = b.ToString( );
          for ( int j = 0; j <= b.isize( ) - n; j++ )
          {    String s = bs.substr( j, n );
               Bool bad = False;
               for ( int k = 0; k < n; k++ )
                    if ( q[j+k] < minq ) bad = True;
               if (bad) continue;
               int p = BinPosition( units, s );
               if ( p >= 0 ) counts[p]++;    }    }
     ReverseSortSync( counts, units );
     int last;
     for ( last = 0; last < counts.isize( ); last++ )
          if ( counts[last] == 0 ) break;
     units.resize(last);

     for ( int pid = 0; pid < (int) pairs.nPairs( ); pid++ )
     {    vec<int> id;
          id.push_back( pairs.ID1(pid), pairs.ID2(pid) );
          vec<String> seq;
          vec<int> passes;
          for ( int z = 0; z < 2; z++ )
          {    int i = id[z];
               for ( int pass = 1; pass <= 2; pass++ )
               {    basevector b = bases[i];
                    if ( pass == 2 ) b.ReverseComplement( );
                    String bs = b.ToString( );
                    for ( int j = 0; j < b.isize( ); j++ )
                    {    for ( int u1 = 0; u1 < units.isize( ); u1++ )
                         {    if ( !bs.Contains( units[u1], j ) ) continue;
                              String s;
                              s = Pic( u1, units.size( ) );
                              while( j <= b.isize( ) - 2*n )
                              {    String x = bs.substr( j + n, n );
                                   int p = Position( units, x );
                                   if ( p < 0 ) break;
                                   s += "." + Pic( p, units.size( ) );
                                   j += n;    }
                              seq.push_back(s);
                              passes.push_back(pass);
                              goto next;     }    }    }
               next: continue;    }
          if ( seq.size( ) == 2 )
          {    if ( passes[0] == 2 ) swap( seq[0], seq[1] );
               cout << seq[0] << " ... " << seq[1] << endl;    }    }    }
