///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// require <cstddef> to use gmp in GCC-4.9
#include <cstddef>
#include <gmpxx.h>

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
// MakeDepend: library GMPXX

#include "CoreTools.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/Lines.h"
#include "paths/long/large/tenx/TenxTools.h"

typedef mpf_class big_float;

void FindTenxNhood( const HyperBasevectorX& hb, const vec<int>& inv,
     const vec<vec<vec<vec<int>>>>& lines, const vec<int>& lens,
     const vec< pair<int,int> >& places, const vec<int>& lhits,
     const vec<double>& bc_frac,
     const vec<int64_t>& all, const vec<uint32_t>& bcs, const int nbc,
     const int nbc_total,
     const vec<String>& genome_names, const vec< vec< pair<int,int> > >& aligns )
{
     // Create index to lines.

     cout << Date( ) << ": making tol4" << endl;
     vec< quad<int,int,int,int> > tol4( hb.E( ) );
     for ( int i = 0; i < lines.isize( ); i++ )
     for ( int j = 0; j < lines[i].isize( ); j++ )
     for ( int k = 0; k < lines[i][j].isize( ); k++ )
     for ( int l = 0; l < lines[i][j][k].isize( ); l++ )
          tol4[ lines[i][j][k][l] ] = make_quad( i, j, k, l );

     // Show lines having most hits.

     {    cout << Date( ) << ": finding lines having most hits" << endl;
          vec< pair<int,int> > lh;
          #pragma omp parallel for
          for ( int i = 0; i < all.isize( ); i++ )
          {    int64_t rid = all[i];
               int eid = places[rid].first;
               if ( eid < 0 ) continue;
               int lid = tol4[eid].first;

               // Ignore if not within 20 kb of line end.

               const vec<vec<vec<int>>>& L = lines[lid];
               int j = tol4[eid].second, k = tol4[eid].third, l = tol4[eid].fourth;
               int left = places[rid].second;
               for ( int m = 0; m < j; m++ )
               {    if ( L[m].empty( ) ) continue;
                    vec<int> lens;
                    for ( int r = 0; r < L[m].isize( ); r++ )
                    {    int n = 0;
                         for ( int u = 0; u < L[m][r].isize( ); u++ )
                              n += hb.Kmers( L[m][r][u] );
                         lens.push_back(n);    }
                    Sort(lens);
                    left += Median(lens);    }
               for ( int u = 0; u < l; u++ )
                    left += hb.Kmers( L[j][k][u] );
               int right = hb.Bases(eid) - places[rid].second - 88;
               for ( int m = j+1; m < L.isize( ); m++ )
               {    if ( L[m].empty( ) ) continue;
                    vec<int> lens;
                    for ( int r = 0; r < L[m].isize( ); r++ )
                    {    int n = 0;
                         for ( int u = 0; u < L[m][r].isize( ); u++ )
                              n += hb.Kmers( L[m][r][u] );
                         lens.push_back(n);    }
                    Sort(lens);
                    right += Median(lens);    }
               for ( int u = l+1; u < L[j][k].isize( ); u++ )
                    right += hb.Kmers( L[j][k][u] );
               if ( left > 20000 && right > 20000 ) continue;

               // Save.

               #pragma omp critical
               {    lh.push( lid, bcs[rid] );
                    lh.push( tol4[ inv[eid] ].first, bcs[rid] );    }    }
          Sort(lh);
          vec< triple<int,int,int> > lhn;
          for ( int i = 0; i < lh.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < lh.isize( ); j++ )
                    if ( lh[j].first != lh[i].first ) break;
               vec<int> b;
               for ( int k = i; k < j; k++ )
                    b.push_back( lh[k].second );
               UniqueSort(b);
               lhn.push( b.size( ), lh[i].first, j - i );
               i = j - 1;    }
          ReverseSort(lhn);
          int prints = 0;
          for ( int i = 0; i < lhn.isize( ); i++ )
          {    int l = lhn[i].second, l2 = -1;
               Bool pair = False;
               if ( i < lhn.isize( ) - 1 )
               {    l2 = lhn[i+1].second;
                    if ( lines[l].front( )[0][0] == inv[ lines[l2].back( )[0][0] ] )
                         pair = True;    }
               vec<String> chrs;
               for ( int j = 0; j < lines[l].isize( ); j++ )
               {    for ( int r = 0; r < lines[l][j].isize( ); r++ )
                    for ( int s = 0; s < lines[l][j][r].isize( ); s++ )
                    {    int e = lines[l][j][r][s];
                         for ( int m = 0; m < aligns[e].isize( ); m++ )
                              chrs.push_back( genome_names[aligns[e][m].first] );
                         e = inv[e];
                         for ( int m = 0; m < aligns[e].isize( ); m++ )
                         {    chrs.push_back( genome_names[aligns[e][m].first] );
                                  }    }    }
               UniqueSort(chrs);
               int obs = lhn[i].first;

               double expect = nbc * bc_frac[l];
               // double expect = double(lhits[l]) * double(nbc)/nbc_total;

               int precision = 100;

               // Let sum = P( X >= obs ) for a a Poisson random variable with
               // lambda = expect.
          
               big_float sum = 0;
               for ( int x = 0; x < obs; x++ )
               {    
                    // Let s1 = expect^x.

                    big_float s1( 1, precision );
                    for ( int j = 0; j < x; j++ )
                         s1 *= expect;

                    // Let s2 = exp(-expect).  Note funny use of 100 term sum.

                    big_float s2( 1, precision );
                    big_float fact( 1, precision );
                    big_float pow( 1, precision );
                    for ( int j = 1; j < 100; j++ )
                    {    pow *= -expect;
                         fact *= j;
                         s2 += pow/fact;    }

                    // Let s = s1 * s2 / x!.

                    big_float s = s1 * s2;
                    for ( int j = 1; j <= x; j++ ) 
                         s /= j;
                    sum += s;    }

               // Let surprise = -log10(1-sum).

               big_float msum = 1 - sum;
               long double surprise = -log10l( msum.get_d( ) );
               if ( lens[l] < 1000 ) continue;
               if ( surprise < 15.0 ) continue;

               cout << "[#" << ++prints << ", " << obs << " barcodes, " 
                    << expect << " expect, "
                    << setprecision(3) << surprise << " surprise"
                    // << ", " << lhn[i].third << " reads"
                    << "] L" << l;
               if ( !pair ) cout << " (l=" << lens[l] << ")";
               else
               {    cout << "/" << l2 << " (l=" << lens[l] << ")";
                    i++;    }
               if ( chrs.nonempty( ) ) cout << " chrs = " << printSeq(chrs);
               cout << "\n";
               if ( prints == 20 ) break;    }    }    }
