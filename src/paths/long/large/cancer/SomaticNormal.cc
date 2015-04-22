///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "math/HoInterval.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/Lines.h"
#include "paths/long/large/cancer/SomaticNormal.h"

void FindBadEdges( const HyperBasevectorX& hb, const vec<int>& inv,
     const vec<vec<vec<vec<int>>>>& lines, const vec< vec< pair<int,int> > >& hits,
     vec<Bool>& bad, const Bool alt )
{    
     vec<int> tol;
     GetTol( hb, lines, tol );
     vec< pair<int,ho_interval> > sn;
     SomaticNormal( sn, alt );
     bad.resize( hb.E( ), False );
     vec<vec<Bool>> flagged( lines.size( ) );
     for ( int l = 0; l < lines.isize( ); l++ )
          flagged[l].resize( lines[l].size( ), False );
     for ( int e = 0; e < hb.E( ); e++ )
     for ( int j = 0; j < hits[e].isize( ); j++ )
     {    int g = hits[e][j].first;
          int start = hits[e][j].second;
          int stop = start + hb.Kmers(e);
          for ( int x = 0; x < sn.isize( ); x++ )
          {    if ( sn[x].first != g ) continue;
               if ( IntervalOverlap( 
                    start, stop, sn[x].second.Start( ), sn[x].second.Stop( ) ) <= 0 )
               {    continue;    }
               bad[e] = True;    
               bad[ inv[e] ] = True;
               int l = tol[e];
               for ( int m = 0; m < lines[l].isize( ); m += 2 )
                    if ( lines[l][m][0][0] == e ) flagged[l][m] = True;    }    }
     for ( int l = 0; l < lines.isize( ); l++ )
     {    for ( int m = 0; m < lines[l].isize( ); m += 2 )
          {    if ( flagged[l][m] )
               for ( int n = m + 2; n <= Min( m + 8, lines[l].isize( ) ); n++ )
               {    if ( flagged[l][n] )
                    {    for ( int r = m + 1; r <= n; r++ )
                         for ( int i = 0; i < lines[l][r].isize( ); i++ )
                         for ( int j = 0; j < lines[l][r][i].isize( ); j++ )
                         {    bad[ lines[l][r][i][j] ] 
                                   = True;    }    }    }    }    }    }
