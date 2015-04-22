///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "FastIfstream.h"
#include "TokenizeString.h"
#include "paths/long/large/svg/Colors.h"
#include "paths/long/large/svg/Svg.h"
#include "random/Random.h"

void GetFields( const String& line, vec<String>& fields )
{    String c;
     Bool in_quote = False;
     for ( int i = 0; i < line.isize( ); i++ )
     {    if ( line[i] == ' ' && !in_quote )
          {    fields.push_back(c);
               c = "";    }
          else if ( line[i] == '"' )
          {    c.push_back( '"' );
               in_quote = !in_quote;    }
          else c.push_back( line[i] );    }
     fields.push_back(c);    }

Bool AtField( const String& x, const String& y )
{    return ( x.Contains( y + "=\"", 0 ) && x.Contains( "\"", -1 ) );    }

String GetField( const String& x, const String& y )
{    return x.Between( y + "=\"", "\"" );    }

double GetFieldDouble( const String& x, const String& y )
{    return x.Between( y + "=\"", "\"" ).Double( );    }

Bool FetchField( const String& x, const String& y, String& z )
{    if ( !AtField( x, y ) ) return False;
     z = GetField( x, y );
     return True;    }

Bool FetchFieldDouble( const String& x, const String& y, double& z )
{    if ( !AtField( x, y ) ) return False;
     z = GetFieldDouble( x, y );
     return True;    }

svg_master::svg_master( const svg_ellipse& x )
{    svgtype_ = "ellipse";
     color_ = x.Color( );
     pos_ = x.Pos( );
     radius_ = x.Radius( );    }

svg_master::svg_master( const svg_text& x )
{    svgtype_ = "text";
     text_ = x.Text( );
     anchor_ = x.Anchor( );
     pos_ = x.Pos( );
     font_ = x.Font( );
     fontsize_ = x.FontSize( );    }

svg_master::svg_master( const svg_polygon& x )
{    svgtype_ = "polygon";
     coords_ = x.Coords( );
     color_ = x.Color( );    }

svg_master::svg_master( const svg_cpath& x )
{    svgtype_ = "cpath";
     coords_ = x.Coords( );
     color_ = x.Color( );
     stroke_width_ = x.StrokeWidth( );
     stroke_opacity_ = x.StrokeOpacity( );    }

void ParseGblock( vec<String>& gblock, String& id, String& sclass, String& title,
     vec<svg_master>& s )
{    for ( int i = 0; i < gblock.isize( ); i++ )
     {    Bool identified = False;
          if ( gblock[i].Contains( "<g", 0 ) )
          {    identified = True;
               String s = gblock[i].After( "<g" );
               String left = s.Before( ">" ), right = s.After( ">" );
               vec<String> x;
               Tokenize( left, ' ', x );
               for ( int j = 0; j < x.isize( ); j++ )
               {    if ( AtField( x[j], "id" ) )
                         id = x[j].Between( "id=\"", "\"" );
                    if ( AtField( x[j], "class" ) )
                         sclass = x[j].Between( "class=\"", "\"" );    }
               title = right.Between( ">", "</title" );    }
          {    svg_cpath p;
               if ( p.Fetch( gblock[i] ) )
               {    identified = True;
                    s.push_back(p);    }    }
          {    svg_polygon p;
               if ( p.Fetch( gblock[i] ) )
               {    identified = True;
                    s.push_back(p);    }    }
          {    svg_ellipse p;
               if ( p.Fetch( gblock[i] ) )
               {    identified = True;
                    s.push_back(p);    }    }
          {    svg_text p;
               if ( p.Fetch( gblock[i] ) )
               {    identified = True;
                    s.push_back(p);    }    }
          if ( !identified ) 
          {    cout << "\nFAILED TO PARSE.\n" << endl;    
               Scram(1);    }    }    }

void ParseSvgFile( const String& fn, vec<String>& head, vec<svg_group>& groups,
     vec<String>& tail, const Bool print_text )
{
     fast_ifstream in(fn);
     String line;

     // Read head.

     getline( in, line );
     if ( in.fail( ) ) ForceAssert( 0 == 1 );
     while( !line.Contains( "<g ", 0 ) )
     {    head.push_back(line);
          getline( in, line );
          if ( in.fail( ) ) ForceAssert( 0 == 1 );    }
     for ( int j = 0; j < 3; j++ )
     {    head.push_back(line);
          getline( in, line );
          if ( in.fail( ) ) ForceAssert( 0 == 1 );    }

     // Read the rest.

     Bool first = True;
     while(1)
     {    if ( !first )
          {    getline( in, line );
               if ( in.fail( ) ) break;    }
          first = False;

          // Read the tail if that's where we're at.

          if ( line == "</g>" )
          {    tail.push_back(line);
               while(1)
               {    getline( in, line );
                    if ( in.fail( ) ) break;
                    tail.push_back(line);    }

               break;     }

          // Discard comments.  Not correct.

          if ( line.Contains( "<!--", 0 ) || line.Contains( "-->" ) ) continue;

          // Read gblock.

          vec<String> gblock;
          while( line != "</g>" )
          {    gblock.push_back(line);
               getline( in, line );
               if ( in.fail( ) ) ForceAssert( 0 == 1 );    }

          // Ignore cluster gblocks.

          if ( gblock[0].Contains( "class=\"cluster\"" ) ) continue;

          // Parse gblock.

          vec<svg_master> s;
          String id, sclass, title;
          ParseGblock( gblock, id, sclass, title, s );
          svg_group g;
          g.SetMaster(s);
          g.SetId(id);
          g.SetTitle(title);
          g.SetClass(sclass);
          if ( !print_text ) g.SetTexting(False);
          groups.push_back(g);    }    }

Bool svg_ellipse::Fetch( String line )
{    if ( !line.Contains( "<ellipse ", 0 ) ) return False;
     line = line.After( "<ellipse " );
     if ( !line.Contains( "/>", -1 ) ) return False;
     line = line.RevBefore( "/>" );
     vec<String> fields;
     GetFields( line, fields );
     for ( int i = 0; i < fields.isize( ); i++ )
     {    const String& f = fields[i];
          if ( f.Contains( "fill=", 0 ) ) continue;
          else if ( FetchField( f, "stroke", color_ ) );
          else if ( FetchFieldDouble( f, "cx", pos_.first ) );
          else if ( FetchFieldDouble( f, "cy", pos_.second ) );
          else if ( FetchFieldDouble( f, "rx", radius_.first ) );
          else if ( FetchFieldDouble( f, "ry", radius_.second ) );
          else return False;    }
     return True;    }

Bool svg_text::Fetch( String line )
{    if ( !line.Contains( "<text ", 0 ) ) return False;
     line = line.After( "<text " );
     if ( !line.Contains( "</text>", -1 ) ) return False;
     line = line.RevBefore( "</text>" );
     if ( !line.Contains( ">" ) ) return False;
     text_ = line.After( ">" );
     line = line.Before( ">" );
     vec<String> fields;
     GetFields( line, fields );
     for ( int i = 0; i < fields.isize( ); i++ )
     {    const String& f = fields[i];
          if ( FetchField( f, "text-anchor", anchor_ ) );
          else if ( FetchField( f, "font-family", font_ ) );
          else if ( FetchFieldDouble( f, "font-size", fontsize_ ) );
          else if ( FetchFieldDouble( f, "x", pos_.first ) );
          else if ( FetchFieldDouble( f, "y", pos_.second ) );
          else return False;    }
     return True;    }

Bool svg_polygon::Fetch( String line )
{    if ( !line.Contains( "<polygon ", 0 ) ) return False;
     line = line.After( "<polygon " );
     if ( !line.Contains( "/>", -1 ) ) return False;
     line = line.RevBefore( "/>" );
     vec<String> fields;
     GetFields( line, fields );
     for ( int i = 0; i < fields.isize( ); i++ )
     {    String& f = fields[i];
          if ( f.Contains( "fill=", 0 ) ) continue;
          else if ( FetchField( f, "stroke", color_ ) );
          else if ( AtField( f, "points" ) )
          {    f = GetField( f, "points" );
               vec<String> x;
               Tokenize( f, ' ', x );
               for ( int j = 0; j < x.isize( ); j++ )
               {    if ( !x[j].Contains( "," ) ) return False;
                    if ( !x[j].Before( "," ).IsDouble( ) ) return False;
                    if ( !x[j].After( "," ).IsDouble( ) ) return False;
                    coords_.push( x[j].Before( "," ).Double( ),
                         x[j].After( "," ).Double( ) );    }    }
          else return False;    }
     return True;    }

Bool svg_cpath::Fetch( String line )
{    if ( !line.Contains( "<path ", 0 ) ) return False;
     line = line.After( "<path " );
     if ( !line.Contains( "/>", -1 ) ) return False;
     line = line.RevBefore( "/>" );
     vec<String> fields;
     GetFields( line, fields );
     stroke_width_ = 1;   // default
     stroke_opacity_ = 1; // default
     for ( int i = 0; i < fields.isize( ); i++ )
     {    String& f = fields[i];
          if ( f.Contains( "fill=", 0 ) ) continue;
          else if ( FetchField( f, "stroke", color_ ) );
          else if ( FetchFieldDouble( f, "stroke-width", stroke_width_ ) );
          else if ( FetchFieldDouble( f, "stroke-opacity", stroke_opacity_ ) );
          else if ( AtField( f, "d" ) )
          {    f = GetField( f, "d" );
               if ( !f.Contains( "M", 0 ) ) return False;
               f = f.After( "M" );
               if ( !f.Contains( "C" ) ) return False;
               String start = f.Before( "C" );
               f = f.After( "C" );
               if ( !start.Contains( "," ) || !start.Before( "," ).IsDouble( )
                    || !start.After( "," ).IsDouble( ) )
               {    return False;    }
               coords_.push( start.Before( "," ).Double( ),
                    start.After( "," ).Double( ) );
               vec<String> x;
               Tokenize( f, ' ', x );
               for ( int j = 0; j < x.isize( ); j++ )
               {    if ( !x[j].Contains( "," ) ) return False;
                    if ( !x[j].Before( "," ).IsDouble( ) ) return False;
                    if ( !x[j].After( "," ).IsDouble( ) ) return False;
                    coords_.push( x[j].Before( "," ).Double( ),
                         x[j].After( "," ).Double( ) );    }    }
          else return False;    }
     return True;    }

pcb_path Cat( const pcb_path& p1, const pcb_path& p2 )
{    pcb_path p = p1;
     p.pop_back( );
     p.append(p2);
     return p;    }

double CB( const double t, const double x1, const double x2, const double x3,
     const double x4 )
{    return (1-t)*(1-t)*(1-t)*x1 + 3*(1-t)*(1-t)*t*x2 + 3*(1-t)*t*t*x3
          + t*t*t*x4;    }

String svg_ellipse::ToString( ) const
{    ostringstream out;
     out << "<ellipse fill=\"" << color_ << "\" stroke=\"" << color_
          << "\" cx=\"" << pos_.first << "\" cy=\"" << pos_.second 
          << "\" rx=\"" << radius_.first << "\" ry=\"" << radius_.second
          << "\"/>\n";
     return out.str( );    }

String svg_text::ToString( ) const
{    ostringstream out;
     out << "<text text-anchor=\"" << anchor_ 
          << "\" x=\"" << pos_.first << "\" y=\"" << pos_.second 
          << "\" font-family=\"" << font_ << "\" font-size=\""
          << fontsize_ << "\">" << text_ << "</text>\n";
     return out.str( );    }

String svg_polygon::ToString( ) const
{    ostringstream out;
     out << "<polygon fill=\"" << color_ << "\" stroke=\"" << color_
          << "\" points=\"";
     for ( int i = 0; i < coords_.isize( ); i++ )
     {    if ( i > 0 ) out << " ";
          out << coords_[i].first << "," << coords_[i].second;    }
     out << "\"/>\n";
     return out.str( );    }

String svg_cpath::ToString( ) const
{    ostringstream out;
     out << "<path fill=\"" << "none" << "\" stroke=\"" << color_
          << "\" stroke-width=\"" << stroke_width_ 
          << "\" stroke-opacity=\"" << stroke_opacity_ << "\" d=\"M";
     for ( int i = 0; i < coords_.isize( ); i++ )
     {    if ( i == 1 ) out << "C";
          if ( i > 1 ) out << " ";
          out << coords_[i].first << "," << coords_[i].second;    }
     out << "\"/>\n";
     return out.str( );    }

String svg_group::ToString( ) const
{    ostringstream out;
     out << "<g id=\"" << Id( ) << "\" class=\"" << Class( )
          << "\"><title>" << Title( ) << "</title>" << endl;
     for ( int i = 0; i < master_.isize( ); i++ )
     {    const svg_master& m = master_[i];
          if ( m.SvgType( ) == "ellipse" )
          {    svg_ellipse x(m);
               out << x.ToString( );    }
          if ( m.SvgType( ) == "text" && print_text_ )
          {    svg_text x(m);
               out << x.ToString( );    }
          if ( m.SvgType( ) == "polygon" )
          {    svg_polygon x(m);
               out << x.ToString( );    }
          if ( m.SvgType( ) == "cpath" )
          {    svg_cpath x(m);
               out << x.ToString( );    }    }
     out << "</g>\n";
     return out.str( );    }

String SvgColor( const vec<int>& c )
{    return "rgb(" + ToString(c[0]) + "%," + ToString(c[1]) + "%,"
          + ToString(c[2]) + "%)";    }

void GetPalette( const int N, vec<vec<int>>& colors, const vec<vec<int>>& s1,
     const double d1, const double d, const int min_color )
{    int fails = 0;
     srandomx(1234567);
     while( colors.isize( ) < N )
     {    vec<int> D(3);
          D[0] = randomx( ) % 100, D[1] = randomx( ) % 100, D[2] = randomx( ) % 100;
          if ( D[0] + D[1] + D[2] < min_color ) continue;
          Bool close = False;
          for ( int j = 0; j < s1.isize( ); j++ )
          {    if ( ColorDist( D, s1[j] ) < d1 )
               {    close = True;
                    break;    }    }
          if (close) 
          {    fails++;
               if ( fails > 1000 )
               {    cout << "1. Could not find enough colors." << endl;
                    PRINT( colors.size( ) );
                    Scram(1);    }
               continue;    }
          for ( int j = 0; j < colors.isize( ); j++ )
          {    if ( ColorDist( D, colors[j] ) < d )
               {    close = True;
                    break;    }    }
          if (close) 
          {    fails++;
               if ( fails > 1000 )
               {    cout << "2. Could not find enough colors." << endl;
                    PRINT( colors.size( ) );
                    Scram(1);    }
               continue;    }
          fails = 0;
          colors.push_back(D);    }    }
