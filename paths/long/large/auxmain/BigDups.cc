///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// BigDups.  Find large perfect duplications in the human genome reference sequence.
//
// See line labeled "not quite right".
//
// Returns near duplicates for unknown reasons.

#include "Bitvector.h"
#include "Intvector.h"
#include "MainTools.h"
#include "Qualvector.h"
#include "TokenizeString.h"
#include "VecUtilities.h"
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/MakeKmerStuff.h"

String chr( const int g )
{    if ( g < 22 ) return ToString( g + 1 );
     if ( g == 22 ) return "X";
     if ( g == 23 ) return "Y";
     return "?";    }

int main( )
{
     RunTime( );

     const int min_dup = 20000;

     cout << Date( ) << ": loading genome" << endl;
     vecbasevector G( "/wga/scr4/bigrefs/grch38/genome.fastb" );
     int ng = 24;
     G.resize(ng);
     vecbasevector GX(G);
     GX.Append(G);
     vecbitvector amb( "/wga/scr4/bigrefs/grch38/genome.lookup.fastamb" );
     amb.resize(ng);
     vecbitvector ambx(amb);
     ambx.Append(amb);
     for ( int g = 0; g < (int) G.size( ); g++ )
     {    GX[ G.size( ) + g ].ReverseComplement( );
          ambx[ G.size( ) + g ].ReverseMe( );    }

     cout << Date( ) << ": making kmers" << endl;
     const int K = 100;
     vec< triple<kmer<K>,int,int> > kmers_plus;
     MakeKmerLookup1( GX, kmers_plus );

     cout << Date( ) << ": finding dups" << endl;
     for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
     {    int64_t j;
          for ( j = i + 1; j < kmers_plus.jsize( ); j++ )
               if ( kmers_plus[j].first != kmers_plus[i].first ) break;
          for ( int64_t k1 = i; k1 < j; k1++ )
          for ( int64_t k2 = i; k2 < j; k2++ )
          {    int g1 = kmers_plus[k1].second, g2 = kmers_plus[k2].second;
               if ( k2 == k1 || !( g1 <= g2 ) || g1 >= (int) G.size( ) ) continue;
               int pos1 = kmers_plus[k1].third, pos2 = kmers_plus[k2].third;
               if ( pos1 > 0 && pos2 > 0 && !ambx[g1][pos1-1]
                    && !ambx[g2][pos2-1] && GX[g1][pos1-1] == GX[g2][pos2-1] )
               {    continue;    }

               if ( g2 == g1 && pos1 > pos2 ) continue;

               // not quite right:
               if ( g2 == g1 + (int) G.size( ) && pos1 > GX[g2].isize( ) - pos2 ) 
                    continue;

               int len = 0;
               for ( int p1 = pos1 + K; p1 < GX[g1].isize( ); p1++ )
               {    len = p1 - pos1;
                    int p2 = pos2 + len;
                    if ( p2 >= GX[g2].isize( ) || GX[g1][p1] != GX[g2][p2]
                         || ambx[g1][p1] || ambx[g2][p2] )
                    {    break;    }    }
               if ( len < min_dup ) continue;
               Bool bad = False;
               for ( int l = 0; l < K; l++ )
               {    if ( ambx[g1][pos1+l] || ambx[g2][pos2+1] ) 
                    {    bad = True;
                         break;    }    }
               if (bad) continue;
               cout << "\n+" << chr(g1) << ":" << pos1 << "-" << pos1 + len
                    << " [len=" << len << "]" << endl;
               if ( g2 < (int) G.size( ) )
                    cout << "+" << chr(g2) << ":" << pos2 << "-" << pos2 + len;
               else
               {    g2 -= G.size( );
                    pos2 = G[g2].isize( ) - pos2 - len;
                    cout << "-" << chr(g2) << ":" << pos2 << "-" << pos2 + len;    }
               cout << endl;    }
          i = j - 1;    }
     cout << "\npeak mem = " << ToStringAddCommas( PeakMemUsageBytes( ) )
          << endl << endl;    }
