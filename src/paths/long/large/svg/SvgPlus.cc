///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "paths/long/large/svg/Svg.h"
#include "paths/long/large/svg/SvgPlus.h"

pcb_path EdgePath( const vec<svg_group_plus>& X, const int e )
{    int v = -1, w = -1;
     pcb_path path;
     for ( int j = 0; j < X.isize( ); j++ )
     {    svg_group_plus g = X[j];
          if ( g.EdgeId( ) == e )
          {    v = g.SourceId( ), w = g.TargetId( );
               break;    }    }
     for ( int j = 0; j < X.isize( ); j++ )
     {    svg_group_plus g = X[j];
          if ( g.VertexId( ) == v )
          {    for ( int k = 0; k < g.Master( ).isize( ); k++ )
               {    svg_master m = g.Master( )[k];
                    if ( m.SvgType( ) == "ellipse" )
                    {    path.push_back( m.Pos( ) );
                         goto found1;    }    }    }    }
     found1:
     for ( int j = 0; j < X.isize( ); j++ )
     {    svg_group_plus g = X[j];
          if ( g.EdgeId( ) == e )
          {    for ( int k = 0; k < g.Master( ).isize( ); k++ )
               {    svg_master m = g.Master( )[k];
                    if ( m.SvgType( ) == "cpath" )
                    {    for ( int l = 1; l < m.Coords( ).isize( ) - 1; l++ )
                              path.push_back( m.Coords( )[l] );
                         goto found2;    }    }    }    }
     found2:
     for ( int j = 0; j < X.isize( ); j++ )
     {    svg_group_plus g = X[j];
          if ( g.VertexId( ) == w )
          {    for ( int k = 0; k < g.Master( ).isize( ); k++ )
               {    svg_master m = g.Master( )[k];
                    if ( m.SvgType( ) == "ellipse" )
                    {    path.push_back( m.Pos( ) );
                         return path;    }    }    }    }
     ForceAssert( 0 == 1 );
     return path;    }

void PointLoc( const vec<svg_group_plus>& X, const int e, const double d,
     double& x, double& y )
{    pcb_path path = EdgePath( X, e );
     double xm = path[0].first + d * ( path.back( ).first - path[0].first );
     for ( int j = 0; j < path.isize( ); j += 3 )
     {    if ( xm >= path[j].first && xm <= path[j+3].first )
          {    for ( double t = 0; t <= 1.0; t += 0.001 )
               {    x = CB( t, path[j].first, path[j+1].first, 
                         path[j+2].first, path[j+3].first );
                    if ( Abs( x - xm ) > 1.0 ) continue;
                    y = CB( t, path[j].second, path[j+1].second, 
                         path[j+2].second, path[j+3].second );
                    return;    }    }    }
     ForceAssert( 0 == 1 );    }
