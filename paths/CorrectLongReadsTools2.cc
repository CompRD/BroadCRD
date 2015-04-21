///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "MainTools.h"
#include "graph/Digraph.h"
#include "paths/BigMapTools.h"
#include "paths/CorrectLongReadsTools2.h"
#include "paths/Useq.h"

void Assess( const vec< vec<int> >& reads, const int nr, const vecbasevector& genome,
     const vecbasevector& genome2, const vecbasevector& unibases, const int K, 
     const int LG, const VecIntPairVec& Glocs,
     const vec<int>& uperfect )
{    cout << "\n" << Date( ) << ": finding perfect reads" << endl;
     vec<Bool> lperfect(nr);
     vec<placementy> places_all;
     vec< vec<ho_interval> > cov( genome.size( ) );
     int nz = 0;
     cout << "\n";
     for ( int id = 0; id < nr; id++ )
     {    const vec<int>& x = reads[id];
          if ( x.empty( ) ) continue;
          nz++;
          basevector b = unibases[ x[0] ];
          for ( int j = 1; j < x.isize( ); j++ )
          {    b.resize( b.isize( ) - (K-1) );
               b = Cat( b, unibases[ x[j] ] );    }
          vec<placementy> places 
               = FindGenomicPlacementsY( 0, b, LG, genome2, Glocs );
          for ( int m = 0; m < places.isize( ); m++ )
          {    placementy p = places[m];
               cov[p.g].push( p.pos, p.Pos );
               int ng = genome2[p.g].isize( ) / 2;
               if ( p.pos >= ng ) continue;
               /*
               cout << m << " placed: " << p.g << "." << p.pos << "-" 
                    << p.Pos << " " << ( p.fw ? "fw" : "rc" ) << "\n";    
               */
                    }
          places_all.append(places);
          static int count(0);
          int pc = 1000;
          if ( places.empty( ) && ++count % pc == 0 )
          {    cout << "[" << count/pc << "] " << id << " imperfect";
               for ( int j = 0; j < x.isize( ); j++ )
               {    cout << " " << x[j];
                    if ( uperfect[ x[j] ] == 0 ) cout << "[FALSE]";    }
               cout << "\n";    }
          lperfect[id] = places.nonempty( );    }
     cout << "\n";
     PRINT(nz);
     cout << PERCENT_RATIO( 3, Sum(lperfect), nz ) << " of reads are "
          << "perfect" << endl;
     vec< vec<ho_interval> > cov2( genome.size( ) ), cov3( genome.size( ) );
     for ( int g = 0; g < (int) genome.size( ); g++ )
     {    Sort( cov[g] );
          for ( int j = 0; j < cov[g].isize( ); j++ )
          {    int k;
               for ( k = j + 1; k < cov[g].isize( ); k++ )
                    if ( cov[g][k].Start( ) > cov[g][j].Start( ) ) break;
               int stop = 0;
               for ( int l = j; l < k; l++ )
                    stop = Max( stop, cov[g][l].Stop( ) );
               ho_interval h( cov[g][j].Start( ), stop );
               cov2[g].push_back(h);
               j = k - 1;    }    }
     for ( int g = 0; g < (int) genome.size( ); g++ )
     {    for ( int j = 0; j < cov2[g].isize( ); j++ )
          {    int k;
               for ( k = j + 1; k < cov2[g].isize( ); k++ )
                    if ( cov2[g][k].Stop( ) > cov2[g][j].Stop( ) ) break;
               cov3[g].push_back( cov2[g][j] );
               j = k - 1;    }    }
     vec<int> overlaps;
     Bool verbose = False;
     for ( int g = 0; g < (int) genome.size( ); g++ )
     {    if (verbose) cout << "\ncoverage of reference " << g << endl;
          for ( int j = 0; j < cov3[g].isize( ); j++ )
          {    if (verbose) cout << cov3[g][j];
               if ( j > 0 ) 
               {    int overlap = cov3[g][j-1].Stop( ) - cov3[g][j].Start( );
                    if (verbose) cout << " [" << overlap << "]";
                    if ( overlap < 1000 ) overlaps.push_back(overlap);    }
               if (verbose) cout << "\n";    
               if ( cov3[g][j].Stop( ) > genome[g].isize( ) + 1000 ) break;    }    }
     int lt0 = 0, lt200 = 0, lt500 = 0, lt640 = 0, lt1000 = 0;
     for ( int i = 0; i < overlaps.isize( ); i++ )
     {    if ( overlaps[i] < 0 ) lt0++;
          else if ( overlaps[i] < 200 ) lt200++;
          else if ( overlaps[i] < 500 ) lt500++;
          else if ( overlaps[i] < 640 ) lt640++;
          else lt1000++;    }
     PRINT5( lt0, lt200, lt500, lt640, lt1000 );    }

void CleanBubbles( const int K, const vecbasevector& unibases, 
     const vec< vec<int> >& nexts, vec< vec< vec<int> > >& upaths, 
     vec< vec<int> >& cores, vec< vec<int> >& ucores, 
     vec< vec< pair<int,int> > >& index, const vecbasevector& genome2, const int LG, 
     const VecIntPairVec& Glocs, const vec< digraphVE<int,int> >& H,
     const vec< digraphVE<int,int> >& Hrc, const vec<int>& uperfect,
     const Bool VERBOSE )
{
     int nuni = unibases.size( ), nr = cores.size( ) / 2;
     cout << Date( ) << ": locating branches" << endl;
     vec< pair< vec<int>, vec<int> > > edits;
     for ( int u = 0; u < nuni; u++ )
     {    if ( nexts[u].size( ) < 2 ) continue;
          vec< vec<int> > exts;
          for ( int id = 0; id < 2*nr; id++ )
          {    for ( int j = 0; j < upaths[id].isize( ); j++ )
               {    for ( int l = 0; l < upaths[id][j].isize( ); l++ )
                    {    if ( upaths[id][j][l] == u ) 
                         {    vec<int> x;
                              for ( int m = l; m < upaths[id][j].isize( ); m++ )
                                   x.push_back( upaths[id][j][m] );
                              exts.push_back(x);    }    }    }    }
          vec< vec<int> > sexts;
          for ( int j = 0; j < exts.isize( ); j++ )
          {    for ( int len = 3; len < exts[j].isize( ); len++ )
               {    vec<int> x = exts[j];
                    x.resize(len);
                    sexts.push_back(x);    }    }
          UniqueSort(sexts);
          for ( int j = 0; j < sexts.isize( ); j++ )
               sexts[j].ReverseMe( );
          Sort(sexts);
          for ( int j = 0; j < sexts.isize( ); j++ )
               sexts[j].ReverseMe( );
          for ( int j = 0; j < sexts.isize( ); j++ )
          {    int k;
               for ( k = j + 1; k < sexts.isize( ); k++ )
                    if ( sexts[k].back( ) != sexts[j].back( ) ) break;
               if ( k - j != 2 )
               {    j = k - 1;
                    continue;    }
               if (sexts[j][1] == sexts[j+1][1] )
               {    j = k - 1;
                    continue;    }
               vec<int> v1 = sexts[j], v2 = sexts[j+1];
               Sort(v1), Sort(v2);
               vec<int> vint = Intersection( v1, v2 );
               if ( vint.size( ) != 2 )
               {    j = k - 1;
                    continue;    }
               int delta = 0;
               for ( int l = j; l < k; l++ )
               {    for ( int m = 0; m < sexts[l].isize( ); m++ )
                    {    int n = unibases[ sexts[l][m] ].isize( ) - K + 1;
                         if ( l == j ) delta += n;
                         else delta -= n;    }    }
               delta = Abs(delta);
               const int max_delta = 50;
               if ( delta > max_delta )
               {    j = k - 1;
                    continue;    }
               vec< vec<int> > v(2);
               v[0] = sexts[j], v[1] = sexts[j+1];
               const int infinity = 1000000000;
               vec<int> score;

               for ( int z = 0; z < index[u].isize( ); z++ )
               {    int id = index[u][z].first; 
                    int p = index[u][z].second, t;
                    if ( ucores[id].Contains( sexts[j], p ) ) t = 0;
                    else if ( ucores[id].Contains( sexts[j+1], p ) ) t = 1;
                    else continue;
                    vec<int> x = cores[id];
                    const digraphVE<int,int>& h = ( id < nr ? H[id] : Hrc[id-nr] );
                    x.SetToSubOf( x, p, v[t].size( ) );
                    int errs = 0, errs_alt = 0;
                    for ( int l = 0; l < x.isize( ) - 1; l++ )
                    {    vec<int> m;
                         for ( int r = 0; r < h.From( x[l] ).isize( ); r++ )
                              if ( h.From( x[l] )[r] == x[l+1] ) m.push_back(r);
                         if ( !m.solo( ) ) // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                         {    PRINT5( id, nr, x[l], x[l+1], m.size( ) ); // XXXXXXXX
                              PRINT2( v[t][l], v[t][l+1] );    } // XXXXXXXXXXXXXXXX
                         ForceAssert( m.solo( ) );
                         errs += h.EdgeObjectByIndexFrom( x[l], m[0] );    }

                    vec< pair< vec<int>, int > > paths;
                    vec<int> y(x);
                    y.resize(1);
                    paths.push( y, 0 );
                    errs_alt = infinity;
                    while( paths.nonempty( ) )
                    {    vec<int> y = paths.back( ).first;
                         int e = paths.back( ).second;
                         paths.pop_back( );
                         if ( y.size( ) == v[1-t].size( ) 
                              && h.Vert( y.back( ) ) == v[1-t].back( ) )
                         {    errs_alt = Min( errs_alt, e );
                              continue;    }
                         int r = y.isize( ) - 1;
                         for ( int m = 0; m < h.From( y[r] ).isize( ); m++ )
                         {    if ( h.Vert( h.From( y[r] )[m] ) == v[1-t][r+1] )
                              {    vec<int> y2(y);
                                   y2.push_back( h.From( y[r] )[m] );
                                   int e2 = e + h.EdgeObjectByIndexFrom( y[r], m );
                                   paths.push( y2, e2 );    }    }    }

                    /*
                    x.resize(1);
                    for ( int r = 0; r < v[1-t].isize( ) - 1; r++ )
                    {    for ( int m = 0; m < h.From( x[r] ).isize( ); m++ )
                         {    if ( h.Vert( h.From( x[r] )[m] ) == v[1-t][r+1] )
                              {    x.push_back( h.From( x[r] )[m] );
                                   errs_alt += h.EdgeObjectByIndexFrom( x[r], m );
                                   break;    }    }
                         if ( x.isize( ) != r + 2 ) break;    }
                    if ( x.size( ) != v[1-t].size( ) ) errs_alt = infinity;
                    */

                    int errs_delta;
                    if ( errs_alt < infinity ) errs_delta = errs_alt - errs;
                    else errs_delta = infinity;
                    if ( t == 0 ) score.push_back(errs_delta);
                    else score.push_back(-errs_delta);    }

               ReverseSort(score);
               int plus = 0, minus = 0;
               int64_t plus2 = 0, minus2 = 0;
               for ( int z = 0; z < score.isize( ); z++ )
               {    if ( score[z] > 0 ) 
                    {    plus++;
                         plus2 += score[z];    }
                    if ( score[z] < 0 ) 
                    {    minus++;
                         minus2 -= score[z];    }    }
               double ratio = double(Max(plus,minus)) / double(Min(plus,minus));
               const double min_ratio = 12.0;
               if (VERBOSE)
               {    cout << "\nbubble group " << edits.size( ) << ":\n";
                    for ( int l = j; l < k; l++ )
                    {    cout << "[" << l-j+1 << "]";
                         for ( int m = 0; m < sexts[l].isize( ); m++ )
                         {    cout << " " << sexts[l][m];
                              if ( genome2.size( ) > 0 
                                   && uperfect[ sexts[l][m] ] == 0 ) 
                              {    cout << "[false]";    }    }
                         basevector b = unibases[ sexts[l][0] ];
                         for ( int m = 1; m < sexts[l].isize( ); m++ )
                         {    b.resize( b.isize( ) - (K-1) );
                              b = Cat( b, unibases[ sexts[l][m] ] );    }
                         if ( genome2.size( ) > 0 )
                         {    vec<placementy> places = FindGenomicPlacementsY( 
                                   0, b, LG, genome2, Glocs );
                              if ( places.empty( ) ) cout << " -- FALSE!";    }
                         cout << "\n";    }    
                    PRINT6( delta, plus, minus, ratio, plus2, minus2 );
                    cout << "winner is " 
                         << ( plus > minus ? "[1]" : "[2]" ) << "\n";    }
               if ( ratio < min_ratio ) continue;
               if ( plus > minus ) edits.push( v[1], v[0] );
               else edits.push( v[0], v[1] );
               /*
               cout << "scores:";
               for ( int z = 0; z < score.isize( ); z++ )
               {    cout << " ";
                    if ( score[z] == infinity ) cout << "inf";
                    else if ( score[z] == -infinity ) cout << "-inf";
                    else cout << score[z];    }
               */
               j = k - 1;    }     }

     if (VERBOSE) cout << "\n";
     for ( int id = 0; id < 2*nr; id++ )
     {    vec<int> &x = ucores[id], &y = cores[id];
          const digraphVE<int,int>& h = ( id < nr ? H[id] : Hrc[id-nr] );
          for ( int p = 0; p < x.isize( ); p++ )
          {    for ( int e = 0; e < edits.isize( ); e++ )
               {    if ( x.Contains( edits[e].first, p ) )
                    {    vec<int> x1, x2, x3, y1, y2, y3;
                         x1.SetToSubOf( x, 0, p );
                         x2 = edits[e].second;
                         x3.SetToSubOf( x, p + edits[e].first.isize( ),
                              x.isize( ) - ( p + edits[e].first.isize( ) ) );
                         y1.SetToSubOf( y, 0, p );
                         y3.SetToSubOf( y, p + edits[e].first.isize( ),
                              y.isize( ) - ( p + edits[e].first.isize( ) ) );

                         vec< vec<int> > Y2x, Y2;
                         vec<int> f;
                         f.push_back( y[p] );
                         Y2x.push_back(f);
                         while( Y2x.nonempty( ) )
                         {    vec<int> a = Y2x.back( );
                              Y2x.pop_back( );
                              int s = a.back( );
                              for ( int j = 0; j < h.From(s).isize( ); j++ )
                              {    if ( h.Vert( h.From(s)[j] ) 
                                        == x[ p + a.isize( ) ] )
                                   {    vec<int> b(a);
                                        b.push_back( h.From(s)[j] );
                                        ForceAssertLe( b.size( ), x2.size( ) );
                                        if ( b.size( ) == x2.size( ) )
                                        {    if ( y3.empty( ) ) Y2.push_back(b);
                                             else if ( h.HasEdge( b.back( ),
                                                  y3.front( ) ) )
                                             {    Y2.push_back(b);    }    }
                                        else Y2x.push_back(b);    }    }    }

                         if (VERBOSE)
                         {    cout << "applying bubble group " << e << " to read "
                                   << id << "\n";    }
                         if ( !Y2.solo( ) )
                         {    if (VERBOSE) PRINT( Y2.size( ) );
                              continue;    }
                         y2 = Y2[0];
                         x = x1, y = y1;
                         x.append(x2), y.append(y2);
                         x.append(x3), y.append(y3);    }    }    }    }
          /*
          y.resize( x.isize( ) );
          for ( int j = 0; j < x.isize( ) - 1; j++ )
          {    vec<int> n;
               if ( id < nr )
               {    for ( int l = 0; l < H[id].From( y[j] ).isize( ); l++ )
                    {    int w = H[id].From( y[j] )[l];
                         if ( H[id].Vert(w) == x[j+1] )
                         {    n.push_back(w);
                              // break;    
                                   }    }    }
               else
               {    for ( int l = 0; l < Hrc[id-nr].From( y[j] ).isize( ); l++ )
                    {    int w = Hrc[id-nr].From( y[j] )[l];
                         if ( Hrc[id-nr].Vert(w) == x[j+1] )
                         {    n.push_back(w);
                              // break;    
                                   }    }    }
               if ( !n.solo( ) )  // ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
               {    PRINT4( j, id, nr, n.size( ) ); // ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
                    cout << "x ="; // ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
                    for ( int l = 0; l < x.isize( ); l++ ) // ZZZZZZZZZZZZZZZZZZZZZZ
                         cout << " " << x[l]; // ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
                    cout << "\n";    } // ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
               ForceAssert( n.solo( ) );
               y[j+1] = n[0];    }    }
          */
     index.clear_and_resize(nuni);
     for ( int id = 0; id < ucores.isize( ); id++ )
     {    const vec<int>& x = ucores[id];
          for ( int j = 0; j < x.isize( ); j++ )
               index[ x[j] ].push( id, j );    }    }

void ExtendAlignmentRight( const vec<int>& x1, const vec<int>& x2,
     const vec<int>& nkmers, const int max_delta, vec< pair<int,int> >& a )
{    int p1 = a.back( ).first, p2 = a.back( ).second;
     if ( p1 == x1.isize( ) - 1 || p2 == x2.isize( ) - 1 ) return;
     int n1 = 0, n2 = 0;
     for ( int j2 = p2 + 1; j2 < x2.isize( ); j2++ )
     {    for ( int j1 = p1 + 1; j1 < x1.isize( ); j1++ )
          {    if ( x1[j1] == x2[j2] )
               {    if ( Abs(n1-n2) <= max_delta )
                    {    a.push( j1, j2 );
                         ExtendAlignmentRight( 
                              x1, x2, nkmers, max_delta, a );    }    }
               n1 += nkmers[ x1[j1] ];    }
          n2 += nkmers[ x2[j2] ];    }    }

void ExtendAlignmentLeft( const vec<int>& x1, const vec<int>& x2,
     const vec<int>& nkmers, const int max_delta, vec< pair<int,int> >& a )
{    int p1 = a.back( ).first, p2 = a.back( ).second;
     if ( p1 == 0 || p2 == 0 ) return;
     int n1 = 0, n2 = 0;
     for ( int j2 = p2 - 1; j2 >= 0; j2-- )
     {    for ( int j1 = p1 - 1; j1 >= 0; j1-- )
          {    if ( x1[j1] == x2[j2] )
               {    if ( Abs(n1-n2) <= max_delta )
                    {    a.push_back( make_pair( j1, j2 ) );
                         ExtendAlignmentLeft( 
                              x1, x2, nkmers, max_delta, a );    }    }
               n1 += nkmers[ x1[j1] ];    }
          n2 += nkmers[ x2[j2] ];    }    }

void ExtendAlignment( const vec<int>& x1, const vec<int>& x2,
     const vec<int>& nkmers, const int p1, const int p2, const int max_delta,
     vec< pair<int,int> >& a )
{    a.clear( );
     a.push( p1, p2 );
     ExtendAlignmentRight( x1, x2, nkmers, max_delta, a );
     a.ReverseMe( );
     ExtendAlignmentLeft( x1, x2, nkmers, max_delta, a );
     a.ReverseMe( );    }

void GetAligns( const int id1, const vec<int>& x1, 
     const vec<int>& nkmers, const vec< vec<int> >& ucores,
     const vec< vec< pair<int,int> > >& cindex,
     vec< pair<ualign, int> >& aligns_ids )
{
     for ( int p1 = 0; p1 < x1.isize( ); p1++ )
     for ( int l = 0; l < cindex[ x1[p1] ].isize( ); l++ )
     {    int id2 = cindex[ x1[p1] ][l].first; 
          if ( id2 == id1 ) continue;
          int p2 = cindex[ x1[p1] ][l].second;
          const vec<int>& x2 = ucores[id2];
          const int max_delta = 10;
          const int min_align = 600;
          useq u1(x1), u2(x2);
          ualign a( u1, u2, p1, p2, max_delta );
          int n = 0;

          for ( int j = 0; j < a.Ties( ).isize( ); j++ )
               n += nkmers[ x1[ a.Tie(j).first ] ];
          if ( a.Ties( ).size( ) < x1.size( ) && n < min_align ) continue;

          if ( a.Ties( ).back( ).first < x1.isize( ) - 1
               && a.Ties( ).back( ).second < x2.isize( ) - 1 )
          {    continue;    }
          if ( a.Ties( ).front( ).first > 0 && a.Ties( ).front( ).second > 0 )
               continue;
          aligns_ids.push( a, id2 );    }
     UniqueSort(aligns_ids);    }

void GetAligns( const int id1, const int p1, const vec<int>& x1, 
     const vec<int>& nkmers, const vec< vec<int> >& ucores,
     const vec< vec< pair<int,int> > >& cindex,
     vec< pair<ualign, int> >& aligns_ids )
{
     for ( int l = 0; l < cindex[ x1[p1] ].isize( ); l++ )
     {    int id2 = cindex[ x1[p1] ][l].first; 
          if ( id2 == id1 ) continue;
          int p2 = cindex[ x1[p1] ][l].second;
          const vec<int>& x2 = ucores[id2];
          const int max_delta = 10;
          const int min_align = 600;
          useq u1(x1), u2(x2);
          ualign a( u1, u2, p1, p2, max_delta );
          int n = 0;

          for ( int j = 0; j < a.Ties( ).isize( ); j++ )
               n += nkmers[ x1[ a.Tie(j).first ] ];
          if ( a.Ties( ).size( ) < x1.size( ) && n < min_align ) continue;

          if ( a.Ties( ).back( ).first < x1.isize( ) - 1
               && a.Ties( ).back( ).second < x2.isize( ) - 1 )
          {    continue;    }
          if ( a.Ties( ).front( ).first > 0 && a.Ties( ).front( ).second > 0 )
               continue;
          aligns_ids.push( a, id2 );    }
     UniqueSort(aligns_ids);    }
