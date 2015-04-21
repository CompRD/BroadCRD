///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// For SVG stuff that is intertangled with assembly stuff.

#ifndef SVG_PLUS_H
#define SVG_PLUS_H

#include "CoreTools.h"
#include "paths/long/large/svg/Svg.h"

class svg_group_plus : public svg_group {

     public:

     svg_group_plus( ) { }
     svg_group_plus( const svg_group& g ) : svg_group(g) { }

     int EdgeId( ) const
     {    for ( int i = 0; i < master_.isize( ); i++ )
          {    if ( master_[i].Text( ) != "" )
               {     if ( master_[i].Text( ).Contains( " " ) ) 
                         return master_[i].Text( ).Before( " " ).Int( );
                     else return master_[i].Text( ).Int( );     }    }
          return -1;    }
     int VertexId( ) const
     {    if ( title_.IsInt( ) ) return title_.Int( );
          return -1;    }
     int SourceId( ) const
     {    String sep = "&#45;&gt;";
          if ( title_.Contains(sep) ) return title_.Before(sep).Int( );
          return -1;    }
     int TargetId( ) const
     {    String sep = "&#45;&gt;";
          if ( title_.Contains(sep) ) return title_.After(sep).Int( );
          return -1;    }

};

// EdgePath: return the 'corrected' poly-cubic-bezier path associated to an edge.
// The correction is that the endpoints of the path are replaced by the centers of
// the vertices at its ends.

pcb_path EdgePath( const vec<svg_group_plus>& X, const int e );

// PointLoc: return the coordinates of a point that is some relative distance d
// (0 <= d <= 1) along a given edge e.  Idiotic implementation and possibly
// completely wrong.

void PointLoc( const vec<svg_group_plus>& X, const int e, const double d,
     double& x, double& y );

#endif
