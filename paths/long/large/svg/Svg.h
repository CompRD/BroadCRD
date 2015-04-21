///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Horrible code for parsing and working with very special kinds of SVG files,
// namely those created from DOT files from DISCOVAR de novo assemblies.  Broken
// in many ways.

#ifndef SVG_H
#define SVG_H

#include "CoreTools.h"
#include "TokenizeString.h"

void GetFields( const String& line, vec<String>& fields );
Bool AtField( const String& x, const String& y );
String GetField( const String& x, const String& y );
double GetFieldDouble( const String& x, const String& y );
Bool FetchField( const String& x, const String& y, String& z );
Bool FetchFieldDouble( const String& x, const String& y, double& z );

class svg_ellipse;
class svg_text;
class svg_polygon;
class svg_cpath;

class svg_master {

     public:

     svg_master( ) { }

     svg_master( const svg_ellipse& x );
     svg_master( const svg_text& x );
     svg_master( const svg_polygon& x );
     svg_master( const svg_cpath& x );

     String SvgType( ) const { return svgtype_; }
     String Color( ) const { return color_; }
     double StrokeWidth( ) const { return stroke_width_; }
     double StrokeOpacity( ) const { return stroke_opacity_; }
     pair<double,double> Pos( ) const { return pos_; }
     pair<double,double> Radius( ) const { return radius_; }
     String Text( ) const { return text_; }
     String Anchor( ) const { return anchor_; }
     String Font( ) const { return font_; }
     double FontSize( ) const { return fontsize_; }
     vec< pair<double,double> > Coords( ) const { return coords_; }

     private:

     String svgtype_;
     String color_;
     double stroke_width_;
     double stroke_opacity_;
     pair<double,double> pos_;
     pair<double,double> radius_;
     String text_;
     String anchor_;
     String font_;
     double fontsize_;
     vec< pair<double,double> > coords_;

};

class svg_ellipse {

     public:

     svg_ellipse( ) { }

     svg_ellipse( const svg_master& x )
     {    color_ = x.Color( );
          pos_ = x.Pos( );
          radius_ = x.Radius( );    }

     Bool Fetch( String line );

     String Color( ) const { return color_; }
     void SetColor( const String& c ) { color_ = c; }
     pair<double,double> Pos( ) const { return pos_; }
     void SetPos( const double x, const double y ) 
     {    pos_.first = x, pos_.second = y; }
     pair<double,double> Radius( ) const { return radius_; }
     void SetRadius( const double rx, const double ry )
     {    radius_.first = rx, radius_.second = ry;    }
     void SetRadius( const double r )
     {    radius_.first = r, radius_.second = r;    }

     String ToString( ) const;

     private:

     String color_;
     pair<double,double> pos_;
     pair<double,double> radius_;

};

class svg_text {

     public:

     svg_text( ) { }

     svg_text( const svg_master& x )
     {    text_ = x.Text( );
          pos_ = x.Pos( );
          anchor_ = x.Anchor( );
          font_ = x.Font( );
          fontsize_ = x.FontSize( );    }

     Bool Fetch( String line );

     String Text( ) const { return text_; }
     pair<double,double> Pos( ) const { return pos_; }
     String Anchor( ) const { return anchor_; }
     String Font( ) const { return font_; }
     double FontSize( ) const { return fontsize_; }

     String ToString( ) const;

     private:

     String text_;
     String anchor_;
     pair<double,double> pos_;
     String font_;
     double fontsize_;

};

class svg_polygon {

     public:

     svg_polygon( ) { }

     svg_polygon( const svg_master& x )
     {    coords_ = x.Coords( );
          color_ = x.Color( );    }

     Bool Fetch( String line );

     String Color( ) const { return color_; }
     vec< pair<double,double> > Coords( ) const { return coords_; }

     String ToString( ) const;

     private:

     vec< pair<double,double> > coords_;
     String color_;

};

class svg_cpath {

     public:

     svg_cpath( ) { }

     svg_cpath( const svg_master& x )
     {    coords_ = x.Coords( );
          color_ = x.Color( );
          stroke_opacity_ = x.StrokeOpacity( );
          stroke_width_ = x.StrokeWidth( );    }

     Bool Fetch( String line );

     String ToString( ) const;

     String Color( ) const { return color_; }
     void SetColor( const String& c ) { color_ = c; }
     vec< pair<double,double> > Coords( ) const { return coords_; }
     void SetCoords( const vec< pair<double,double> >& c ) { coords_ = c; }
     double StrokeWidth( ) const { return stroke_width_; }
     void SetStrokeWidth( const double w ) { stroke_width_ = w; }
     double StrokeOpacity( ) const { return stroke_opacity_; }
     void SetStrokeOpacity( const double o ) { stroke_opacity_ = o; }

     private:

     vec< pair<double,double> > coords_;
     String color_;
     double stroke_width_;
     double stroke_opacity_;

};

class svg_group {

     public:

     svg_group( ) : print_text_(True) { }

     vec<svg_master> Master( ) const { return master_; }
     String Class( ) const { return class_; } // should be "node" or "edge"
     String Id( ) const { return id_; }
     String Title( ) const { return title_; }
     void SetMaster( const vec<svg_master>& x ) { master_ = x; }
     void SetId( const String& x ) { id_ = x; }
     void SetClass( const String& x ) { class_ = x; }
     void SetTitle( const String& x ) { title_ = x; }
     void SetTexting( const Bool& x ) { print_text_ = x; }

     String ToString( ) const;

     protected:

     vec<svg_master> master_;
     String id_;
     String class_;
     String title_;
     Bool print_text_;

};

void ParseGblock( vec<String>& gblock, String& id, String& sclass, String& title,
     vec<svg_master>& s );

void ParseSvgFile( const String& fn, vec<String>& head, vec<svg_group>& groups,
     vec<String>& tail, const Bool print_text );

// A pcb_path is a poly-cubic-bezier path.  Actually as used, all that matters is
// that it is a list of points p1,...,pn that describes a curve from p1 to pn.

class pcb_path : public vec< pair<double,double> > {
     public:
     pcb_path( ) { }
     pcb_path( const vec< pair<double,double> >& p ) 
          : vec< pair<double,double> >(p) { }
};

pcb_path Cat( const pcb_path& p1, const pcb_path& p2 );

// cubic Bezier

double CB( const double t, const double x1, const double x2, const double x3,
     const double x4 );

String SvgColor( const vec<int>& c );

// Return N colors of distance at least d1 from colors in s1 and distance at
// least d from each other, and no less than min_color (to avoid dark colors).  
// May fail.

void GetPalette( const int N, vec<vec<int>>& colors, const vec<vec<int>>& s1,
     const double d1, const double d, const int min_color );

#endif
