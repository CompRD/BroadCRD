///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "VecUtilities.h"
#include "paths/long/ReadStack.h"
#include "paths/long/large/explore/SnorgleTools.h"

void FindPaths( const readstack& s, vec<basevector>& p, const Bool FP_LOGGING )
{    
     // Heuristics.

     const int K = 20;
     const int max_paths = 100;
     const int min_ratio = 4;
     const int min_qsum = 80;

     p.clear( );
     vec<basevector> paths;

     if (FP_LOGGING)
     {    cout << "\nFP_LOGGING:";
          for ( int i = 0; i < s.Cols( ); i++ )
          {    if ( i % 80 == 0 ) cout << "\n";
               if ( !s.Def( 0, i ) ) cout << "-";
               else cout << as_base( s.Base( 0, i ) );    }
          cout << "\n";
          for ( int i = 0; i < s.Cols( ); i++ )
          {    if ( i % 80 == 0 ) cout << "\n";
               if ( !s.Def( 1, i ) ) cout << "-";
               else cout << as_base( s.Base( 1, i ) );    }
          cout << "\n\n";    }

     if ( s.Cols( ) < K ) return;

     // Set up initial kmers.  If the read agrees with the consensus, use that.
     // Otherwise, find initial kmers that appear at least twice in the reads.  
     // These constitute the initial paths.

     basevector jcon = s.Consensus1( );
     Bool mismatch = False;
     for ( int i = 0; i < K; i++ )
     {    if ( !s.Def(0,i) || s.Base(0,i) != jcon[i] )
          {    mismatch = True;    }    }
     if ( !mismatch )
     {    basevector b;
          for ( int i = 0; i < K; i++ )
               b.push_back( s.Base( 0, i ) );
          paths.push_back(b);    }
     else
     {    vec<basevector> bx;
          for ( int j = 0; j < s.Rows( ); j++ )
          {    Bool def = True;
               for ( int i = 0; i < K; i++ )
               {    if ( !s.Def( j, i ) )
                    {    def = False;
                         break;    }    }
               if ( !def ) continue;
               basevector b;
               for ( int i = 0; i < K; i++ )
                    b.push_back( s.Base( j, i ) );
               bx.push_back(b);    }
          Sort(bx);
          for ( int i = 0; i < bx.isize( ); i++ )
          {    int j = bx.NextDiff(i);
               if ( j - i >= 2 ) paths.push_back( bx[i] );
               i = j - 1;    }    }
     if (FP_LOGGING) cout << "initial paths: " << paths.size( ) << endl;

     // Extend paths.  At each step, a path has to have at least two extensions.
     // The best extension can beat another according to certain qsum criteria.

     for ( int j = K; j < s.Cols( ); j++ )
     {    vec<basevector> paths2;
          for ( int m = 0; m < paths.isize( ); m++ )
          {    const basevector& b = paths[m];
               vec<int> ext( 4, 0 ), qsum( 4, 0 ), ids( 4, vec<int>::IDENTITY );
               vec<int> qsum2( 4, 0 );
               for ( int r = 0; r < s.Rows( ); r++ )
               {    Bool mismatch = False;
                    for ( int i = 0; i < K - 1; i++ )
                    {    if ( s.Base( r, j - K + 1 + i ) != b[ j - K + 1 + i ] )
                         {    mismatch = True;
                              break;    }    }
                    if (mismatch) continue;
                    if ( !s.Def( r, j ) ) continue;
                    ext[ s.Base( r, j ) ]++;
                    qsum[ s.Base( r, j ) ] += s.Qual( r, j );
                    if ( s.Qual(r,j) > 2 ) 
                         qsum2[ s.Base( r, j ) ] += s.Qual( r, j );    }
               // cout << "m = " << m << ", ext = " << printSeq(ext) << endl; // XXX
               ReverseSortSync( qsum, qsum2, ext, ids );
               for ( int l = 0; l < 4; l++ )
               {    if ( ext[l] < 2 ) continue;
                    if ( l > 0 && qsum[0] >= min_qsum 
                         && qsum[0] >= min_ratio * qsum[l] )
                    {    continue;    }
                    if ( l > 0 && qsum2[0] >= min_qsum && qsum2[0] 
                         >= min_ratio * qsum2[l] )
                    {    continue;    }

                    /*
                    if ( l > 0 ) // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                         PRINT6(j, l, qsum[0], qsum[l], qsum2[0], qsum2[l]); // XXXX
                    */

                    basevector b2(b);
                    b2.push_back( ids[l] );
                    paths2.push_back(b2);    }    }
          paths = paths2;
          if ( paths.isize( ) > max_paths ) return;    }

     // Clean paths.

     for ( int i = 0; i < paths.isize( ); i++ )
     {    Bool mismatch = False;
          for ( int j = 0; j < s.Cols( ); j++ )
          {    for ( int m = 0; m < 2; m++ )
               {    if ( s.Qual( m, j ) >= 30 && paths[i][j] != s.Base( m, j ) )
                    {    mismatch = True;
                         break;    }    }
               if ( mismatch ) break;    }
          if ( !mismatch ) p.push_back( paths[i] );    }    }
