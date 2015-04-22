///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// mins.size( ) = 97
// ambs = 15
// median = 8.54785
// max = 81.2471
// sum = 1119.46

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"
#include "paths/long/large/Lines.h"

int main( int argc, char** argv )
{     
     RunTime( );
     BeginCommandArgumentsNoHeader;
     CommandArgument_String(FID);
     EndCommandArguments;

     /*
     fid = 20; // 9, max = 21.6%, one ambiguity
     fid = 21; // 8, max = 14.4%
     fid = 22; // 4, max = 21.1%
     fid = 23; // 2, max = 16.0%
     fid = 25; // 2, max = 24.0%
     fid = 26; // 2, max = 8.1%
     fid = 27; // 1, max = 4.6%
     fid = 28; // 1, max = 14.2%
     fid = 29; // 1, max = 12.5%
     fid = 30; // 2, max = 20.9%, one ambiguous
     fid = 31; // 4, max = 30.9%, one ambiguous
     fid = 32; // 1, max = 22.2%
     fid = 33; // 2, max = 12.4%, one ambiguous
     fid = 34; // 5, max = 39.1%, two ambiguous
     fid = 35; // 4, max = 40.2%, one ambiguous // loop=189 kmers, counting problem
     fid = 36; // 3, max = 12.4%
     fid = 37; // 2, max = 15.5%, one ambiguous
     fid = 38; // 3, max = 20.1%, one ambiguous
     fid = 39; // 3, max = 27.8%, one ambiguous
     */

     vec<int> fids;
     if ( FID == "all" )
     {    for ( int i = 1; i <= 40; i++ )
               fids.push_back(i);    }
     else if ( FID.IsInt( ) ) fids.push_back( FID.Int( ) );

     vec<double> mins;
     int ambs = 0;

     for ( int fi = 0; fi < fids.isize( ); fi++ )
     {
     int fid = fids[fi];
     cout << "\n";
     PRINT(fid);

     String work_dir = "/wga/scr4/jaffe/GapToy/FDA/" + ToString(fid);
     HyperBasevector hb;
     BinaryReader::readFile( work_dir + "/a.final/a.hbv", &hb );
     vec<int> inv;
     BinaryReader::readFile( work_dir + "/a.final/a.inv", &inv );
     ReadPathVec paths( work_dir + "/a.final/a.paths" );

     VecULongVec paths_index;
     invert( paths, paths_index );

     vec<int> to_right;
     hb.ToRight(to_right);

     vec< vec< vec< vec<int> > > > lines;
     BinaryReader::readFile( work_dir + "/a.final/a.lines", &lines );

     vec<int> lens, tol;
     GetLineLengths( hb, lines, lens );
     GetTol( hb, lines, tol );

     // Heuristics.

     const int min_min_flank = 1000;

     // Find inverse lines.

     int count = 0;
     for ( int v = 0; v < hb.N( ); v++ )
     {    if ( !hb.From(v).solo( ) ) continue;
          if ( hb.To(v).size( ) != 2 ) continue;
          int l1 = tol[ hb.IFrom(v,0) ];
          int l1b = lines[l1].back( )[0][0];

          // Skipping rc for now.

          if ( tol[ inv[l1b] ] < l1 ) continue;

          int w = to_right[l1b];
          if ( !hb.To(w).solo( ) ) continue;
          if ( hb.From(w).size( ) != 2 ) continue;
          for ( int j = 0; j < 2; j++ )
          {    int l2 = tol[ hb.IFrom(w,j) ];
               int l2b = lines[l2].back( )[0][0];

               if ( to_right[l2b] != v ) continue;
               int m1, m2 = tol[ hb.IFrom(w,1-j) ];
               if ( hb.ITo(v,0) == l2b ) m1 = tol[ hb.ITo(v,1) ];
               else m1 = tol[ hb.ITo(v,0) ];

               vec<int64_t> pl1 = LinePids( lines[l1], hb, inv, paths, paths_index );
               vec<int64_t> pl2 = LinePids( lines[l2], hb, inv, paths, paths_index );
               vec<int64_t> pm1 = LinePids( lines[m1], hb, inv, paths, paths_index );
               vec<int64_t> pm2 = LinePids( lines[m2], hb, inv, paths, paths_index );

               if ( Max( lens[m1], lens[m2] ) < min_min_flank ) continue;

               cout << "\n#" << ++count << ". ";
               cout << "l1b = " << l1b << "  (l1 = " << l1 << ", l2 = " << l2
                    << ", m1 = " << m1 << ", m2 = " << m2 << ")" << endl;
               PRINT4( pl1.size( ), pl2.size( ), pm1.size( ), pm2.size( ) );

               double cl1 = double(pl1.size()) / lens[l1];
               double cl2 = double(pl2.size()) / lens[l2];
               double cm1 = double(pm1.size()) / lens[m1];
               double cm2 = double(pm2.size()) / lens[m2];

               int mref = ( lens[m1] >= lens[m2] ? m1 : m2 );

               const int min_flank = 10000;
               int ladd = hb.K( ) - 1 - 60;
               // int ladd = 0;
               int e1 = lines[m1].back( )[0][0];
               if ( hb.Kmers(e1) >= min_flank )
               {    vec<int64_t> p = Pids( e1, hb, inv, paths, paths_index );
                    cm1 = double( p.size( ) ) / ( hb.Kmers(e1) + ladd );    }

               int e2 = lines[m2].front( )[0][0];
               if ( hb.Kmers(e2) >= min_flank )
               {    vec<int64_t> p = Pids( e2, hb, inv, paths, paths_index );
                    cm2 = double( p.size( ) ) / ( hb.Kmers(e2) + ladd );    }

               double cmref = ( mref == m1 ? cm1 : cm2 );

               if ( lens[m1] >= min_flank && lens[m2] >= min_flank )
                    cmref = ( cm1 + cm2 ) / 2.0;

               PRINT4( cl1, cl2, cm1, cm2 );
               PRINT3( cmref, lens[l1], lens[l2] );

               vec<int64_t> x = pl1;
               x.append(pl2);
               UniqueSort(x);

               double err, r, c;
               vec<double> errs(7);
               double min_errs = 1000000000;
               for ( int d = 1; d <= 6; d++ )
               {    if ( d == 1 ) r = cl1/cmref;
                    else
                    {    int len1 = lens[l1];
                         int len2 = lens[l2];
                         c = double( x.size( ) ) / ( d*len1 + (d-1)*len2 + ladd );
                         r = c/cmref;    }
                    if ( r >= 1 ) err = 100 * (r-1);
                    else err = 100 * ( 1/r - 1 );
                    errs[d] = err;
                    min_errs = Min( errs[d], min_errs );
                    if ( errs[d] >= 20 && errs[d] > 2 * min_errs ) 
                    {    errs.resize(d+1);
                         break;    }    }
               int prints = 0;
               for ( int d = 1; d < errs.isize( ); d++ )
               {    if ( errs[d] >= 20 && errs[d] > 2 * min_errs ) continue;
                    cout << d << ". " << errs[d] << "%" << endl;
                    prints++;    }
               mins.push_back(min_errs);
               if ( prints > 1 ) ambs++;    }    }    }

     cout << "\n";
     cout << "SUMMARY STATS\n\n";
     PRINT( mins.size( ) );
     PRINT(ambs);
     Sort(mins);
     if ( mins.nonempty( ) )
     {    cout << "median = " << Median(mins) << endl;
          cout << "max = " << Max(mins) << endl;
          cout << "sum = " << Sum(mins) << endl;    }

     }
