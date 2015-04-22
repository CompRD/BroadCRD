///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef SUPER_SCORE_H
#define SUPER_SCORE_H

#include "CoreTools.h"

auto fsuper = []( int x )
{    if ( x <= 400 ) return 2.5;
     else if ( x > 400 && x <= 1000 ) return 1000.0/x;
     else if ( x > 1000 && x <= 2000 ) return 4 - 3*x/1000.0;
     else return -2.0;    };

// Note asymmetric:

auto gsuper = [&]( int i1, int i2, int j1, int j2 )
{    vec< pair<int,int> > r, l;
     int s1 = 0, s2 = 0;
     for ( int k = j1 + 1; k < X[i1].isize( ); k++ )
     {    s1 += X[i1][k];
          r.push( s1, 1 );    }
     for ( int k = j2 + 1; k < X[i2].isize( ); k++ )
     {    s2 += X[i2][k];
          if ( s2 > s1 ) break;
          r.push( s2, 2 );    }
     s1 = 0, s2 = 0;
     for ( int k = j1 - 1; k >= 0; k-- )
     {    s1 += X[i1][k];
          l.push( s1, 1 );    }
     for ( int k = j2 - 1; k >= 0; k-- )
     {    s2 += X[i2][k];
          if ( s2 > s1 ) break;
          l.push( s2, 2 );    }
     Sort(r), Sort(l);
     double score = 0;
     for ( int pass = 1; pass <= 2; pass++ )
     {    const vec< pair<int,int> >& x = ( pass == 1 ? r : l );
          for ( int i = 0; i < x.isize( ); i++ )
          {    if ( x[i].second != 2 ) continue;
               int best = 2000;
               for ( int j = i - 1; j >= 0; j-- )
               {    if ( x[i].first - x[j].first >= best ) break;
                    if ( x[j].second != 1 ) continue;
                    best = x[i].first - x[j].first;    }
               for ( int j = i + 1; j < x.isize( ); j++ )
               {    if ( x[j].first - x[i].first >= best ) break;
                    if ( x[j].second != 1 ) continue;
                    best = x[j].first - x[i].first;    }
               score += fsuper(best);    }    }
     return score;    };

#endif
