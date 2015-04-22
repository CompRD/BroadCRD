///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "VecUtilities.h"
#include "paths/long/large/bng/BngAlign.h"
#include "paths/long/large/bng/BngAlignCore.h"

double RelDiff( const double d1, const double d2 )
{    return Abs( d1 - d2 ) / Max( d1, d2 );    }

void BngAlignCore( const vec<vec<double>>& X,
     const vec< quad<int,int,int,int> >& zmatches,
     const int64_t low, const int64_t high,
     const vec<int>& S,
     const Bool ALIGN_LOGGING, const Bool SHOW_ALL_SEEDS, const Bool ALIGN_DETAILS,
     const Bool ONE_SEED, const Bool RFILTER, int& PRINT_FAIL,
     int& align_calls,
     vec< triple< int, int, vec<int> > >& aligns, vec<double>& scores,
     double& ac1, double& ac2, double& ac3, double& ac4, double& ac5,
     const int min_direct, const double max_score, const double max_total_err,
     const int s1, const int i1 )

{

#include "paths/long/large/bng/SuperScore.h"

     double bclock1 = 0, bclock2 = 0, bclock3 = 0, bclock4 = 0, bclock5 = 0;
     bclock1 = -WallClockTime( );

     // Parallel structures.

     vec< triple< int, int, vec<int> > > alignsi;
     vec<String> areports;
     vec< pair<int,int> > seed_points;
     vec<double> scoresi;

     int align_callsi = 0;
     
     // Set up logging.

     ostringstream xout;
     if (ALIGN_LOGGING)
     {    xout << "\n==================================================="
               << "=================================\n\n";
          xout << "ALIGNMENTS OF READ " << S[s1] << endl;    }
     
     // Find alignments.
     
     xout << "\nzmatches says there are " << high - low << " hits" << endl;
     vec< triple<int,int,int> > hits;
     for ( int64_t j = low; j < high; j++ )
     {    int i2 = zmatches[j].second;
          int j1 = zmatches[j].third, j2 = zmatches[j].fourth;
          hits.push( i2, j1, j2 );    }

     if (ALIGN_LOGGING)
     {    vec<int> zreads;
          for ( int64_t j = low; j < high; j++ )
          {    int64_t k;
               for ( k = j + 1; k < high; k++ )
                    if ( zmatches[k].second != zmatches[j].second ) break;
               zreads.push_back( zmatches[j].second );
               j = k - 1;    }
          xout << "involving " << zreads.size( ) << " reads" << endl;    }
          // xout << "= " << printSeq(zreads) << endl;

     bclock1 += WallClockTime( );
     for ( int h1 = 0; h1 < hits.isize( ); h1++ )
     {    bclock2 -= WallClockTime( );
          int i2 = hits[h1].first;
          int h2;
          for ( h2 = h1 + 1; h2 < hits.isize( ); h2++ )
               if ( hits[h2].first != hits[h1].first ) break;

          const vec<double> &x1 = X[i1], &x2 = X[i2];
          vec<int> y1, y2;
          for ( int j = 0; j < x1.isize( ); j++ )
               y1.push_back( x1[j] );
          for ( int j = 0; j < x2.isize( ); j++ )
               y2.push_back( x2[j] );

          vec< pair<int,int> > seeds; // (j1,j2)
          for ( int h = h1; h < h2; h++ )
               seeds.push( hits[h].second, hits[h].third );

          bclock2 += WallClockTime( );
          bclock3 -= WallClockTime( );

          // Compute offsets.

          vec< triple<int,int,int> > so;
          for ( int u = 0; u < seeds.isize( ); u++ )
          {    int j1 = seeds[u].first, j2 = seeds[u].second;
               int offset = 0;
               for ( int a = 0; a < j1; a++ )
                    offset += y1[a];
               for ( int a = 0; a < j2; a++ )
                    offset -= y2[a];
               so.push( offset, j1, j2 );    }
          Sort(so);

          bclock3 += WallClockTime( );
          bclock4 -= WallClockTime( );

          // Filter seeds more.

          if ( so.nonempty( ) )
          {    const int min_seed2 = 3000;
               const int min_seeds2 = 4;
               const double prox2 = 0.1;
               const int seed_radius2 = 500;
               const int min_span2 = 40000;
               vec< triple<int,int,int> > so2;
               for ( int u = 0; u < so.isize( ); u++ )
               {    int offset = so[u].first; 
                    int j1 = so[u].second, j2 = so[u].third;
                    vec< pair<int,int> > matches;
                    for ( int l1 = 0; l1 < y1.isize( ); l1++ )
                    {    if ( y1[l1] < min_seed2 ) continue;
                         for ( int l2 = 0; l2 < y2.isize( ); l2++ )
                         {    if ( y2[l2] < min_seed2 ) continue;
                              if ( RelDiff( y1[l1], y2[l2] ) > prox2 ) 
                                   continue;
                              int offsetx = 0;
                              for ( int a = 0; a < l1; a++ )
                                   offsetx += y1[a];
                              for ( int a = 0; a < l2; a++ )
                                   offsetx -= y2[a];
                              if ( Abs( offset - offsetx ) 
                                   <= seed_radius2 )
                              {    matches.push( l1, l2 );    }    }    }
                    if ( matches.isize( ) >= min_seeds2 ) 
                    {    int span = 0;
                         for ( int l = matches.front( ).first;
                              l <= matches.back( ).first; l++ )
                         {    span += y1[l];    }
                         if ( span >= min_span2 )
                              so2.push_back( so[u] );    }    }
               so = so2;    }
     
          // Show all seeds.
                         
          if ( SHOW_ALL_SEEDS && so.nonempty( ) )
          {    xout << "\n";
               PRINT2_TO( xout, i1, i2 );
               for ( int u = 0; u < so.isize( ); u++ )
               {    int offset = so[u].first; 
                    int j1 = so[u].second, j2 = so[u].third;
                    PRINT3_TO( xout, offset, j1, j2 );     }    }
          bclock4 += WallClockTime( );

          // Let's just use one seed.

          if ( ONE_SEED && so.nonempty( ) )
          {    vec< triple<int,int,int> > so3;
               so3.push_back( so[ so.size( ) / 2 ] );
               so = so3;    }

          // Start alignment loop.

          for ( int u = 0; u < so.isize( ); u++ )
          {    int j1 = so[u].second, j2 = so[u].third;
               int offset = so[u].first;

               // Align.

               int rstart, rstop;
               int ndirect;
               vec<int> p;
               int sum1, sum2;
               bclock5 -= WallClockTime( );
               ostringstream out;
               double score = BngAlign( y1, y2, j1, j2, out, p, 
                    rstart, rstop, ndirect, sum1, sum2, 2 );
               bclock5 += WallClockTime( );
               align_callsi++;
               int len = rstop - rstart + 1;
               double total_err = RelDiff( sum1, sum2 );

               // Test alignment quality.

               Bool fail = False;
               if ( ndirect < min_direct || score > max_score ) 
                    fail = True;
               if ( total_err > max_total_err ) fail = True;
               if ( fail && PRINT_FAIL > 0 )
               {
                    #pragma omp critical
                    {    PRINT_FAIL--;
                         cout << "\nFAIL\n";
                         PRINT4( rstart, rstop, len, ndirect );
                         cout << i1 << "[l=" << y1.size( ) << "]." << j1
                              << "(" << y1[j1] << ")"
                              << " aligned to " << i2 << "[l=" 
                              << y2.size( ) 
                              << "]." << j2 << "(" << y2[j2] << ")"
                              << ", offset = " << offset << ", score = " 
                              << score << ", total_err = " 
                              << total_err << endl;
                         cout << out.str( );    

                         cout << "Rfull: " << printSeq(y1) << endl;
                         cout << "Afull: " << printSeq(y2) << endl;
                         double super = gsuper( i1, i2, j1, j2 );
                         cout << "super score = " 
                              << super << endl;    }    }
               if (fail) continue;

               // Report alignment.

               alignsi.push( i1, i2, p );
               seed_points.push( j1, j2 );
               scoresi.push_back(score);
               if (ALIGN_LOGGING)
               {    ostringstream rout;
                    PRINT4_TO( rout, rstart, rstop, len, ndirect );
                    rout << i1 << "[l=" << y1.size( ) << "]"
                         << " aligned to " << i2 << "[l=" << y2.size( ) 
                         << "]" << ", offset = " << offset 
                         << ", score = " << score 
                         << ", total_err = " << total_err << endl;    
                    if (ALIGN_DETAILS) rout << out.str( );
                    rout << "super score = " << gsuper( i1, i2, j1, j2 ) 
                         << endl;
                    areports.push_back( rout.str( ) );    }    }

          h1 = h2 - 1;    }

     // Report and save alignments.

     if (ALIGN_LOGGING) 
          SortSync( alignsi, areports, seed_points, scoresi );
     else SortSync( alignsi, seed_points, scoresi );

     vec<double> rscores;
     for ( int i = 0; i < alignsi.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < alignsi.isize( ); j++ )
               if ( alignsi[j] != alignsi[i] ) break;
          rscores.push_back( scoresi[i] * 100 );
          i = j - 1;    }
     Sort(rscores);

     if ( RFILTER && alignsi.nonempty( ) )
     {    const int score_batch = 2;
          const double score_add = 1.5;
          double mscore 
               = rscores[ Min( rscores.isize( ) - 1, score_batch - 1 ) ];
          vec<Bool> pdel( alignsi.isize( ), False );
          for ( int i = 0; i < alignsi.isize( ); i++ )
               if ( 100 * scoresi[i] > mscore + score_add ) pdel[i] = True;
          EraseIf( alignsi, pdel );
          if (ALIGN_LOGGING) EraseIf( areports, pdel );
          EraseIf( seed_points, pdel );
          EraseIf( scoresi, pdel );    }

     rscores.clear( );
     for ( int i = 0; i < alignsi.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < alignsi.isize( ); j++ )
               if ( alignsi[j] != alignsi[i] ) break;
          rscores.push_back( scoresi[i] * 100 );
          i = j - 1;    }
     Sort(rscores);
     
     vec<Bool> adel( alignsi.size( ), False );
     if (ALIGN_LOGGING)
     {    xout << "\nscores(%): " << setprecision(3)
               << printSeq(rscores) << endl;    }
     for ( int i = 0; i < alignsi.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < alignsi.isize( ); j++ )
               if ( alignsi[j] != alignsi[i] ) break;
          if (ALIGN_LOGGING)
          {    xout << "\nseed points:";
               for ( int k = i; k < j; k++ )
               {    int j1 = seed_points[k].first, 
                         j2 = seed_points[k].second;
                    int i1 = alignsi[i].first, i2 = alignsi[i].second;
                    xout << " (" << j1 << "," << j2 << ")" << " = (" 
                         << int(X[i1][j1]) << "," << int(X[i2][j2]) 
                         << ")";    }
               xout << "\n" << areports[i];    }
          for ( int k = i + 1; k < j; k++ )
               adel[k] = True;
          i = j - 1;    }
     EraseIf( alignsi, adel );
     EraseIf( scoresi, adel );
     #pragma omp critical
     {    if (ALIGN_LOGGING) cout << xout.str( );    
          ac1 += bclock1, ac2 += bclock2, ac3 += bclock3, ac4 += bclock4;
          ac5 += bclock5;
          aligns.append(alignsi);    
          align_calls += align_callsi;    
          scores.append(scoresi);    }    }
