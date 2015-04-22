///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Example: this shows the four chromosomes in HCC1143 (tumor only) on top of
// the HCC1143+BL assembly, with uniquely aligned 10X reads shown as dots.
// This will only work so long as the particular HCC1143+BL assembly is alive.
//
// To get this into a powerpoint document, first print to pdf at 40% and without
// headers and footers.  Bring into powerpoint and crop vertically.  Then scale
// to 125% (from 68% which is what it will be).

#include "FastIfstream.h"
#include "MainTools.h"
#include "VecUtilities.h"
#include "math/Functions.h"
#include "paths/long/large/svg/Svg.h"
#include "paths/long/large/svg/SvgPlus.h"
#include "random/Random.h"
#include "random/Shuffle.h"
#include "paths/long/large/tenx/TenxDirs.h"

int main( int argc, char* argv[] )
{
     RunTime( );

     String default_seeds = "10258506,5332247,10269394";
     String default_edges = "5246658,10258506,10254768,5332247,10269394,10265790,5412188";
     int default_dataset = 1;

     BeginCommandArguments;
     CommandArgument_Int_OrDefault_Doc(N, default_dataset,
             "assembly number to use (10X numbering 1-6)");
     CommandArgument_String_OrDefault_Doc(DIR, ".",
             "output and input for the SVG");
     CommandArgument_Bool_OrDefault_Doc(CACHED, True,
             "whether to use cached Nhood/EdgeInfo");
     CommandArgument_String_OrDefault_Doc(SEEDS, default_seeds,
             "seeds for NhoodInfo");
     CommandArgument_String_OrDefault_Doc(EDGES, default_edges,
             "path edges for EdgeInfo10X");

     CommandArgument_Int_Doc(DISPLAY, "series of displays [1,9]");
     EndCommandArguments;


     // Define series of displays.

     int display = DISPLAY;
     // display = 1; // just the bare graph
     // display = 2; // just the rearranged path with one barcode
     // display = 3; // just the rearranged path



     // display = 7; // everything
     // display = 8; // drop the dots
     // display = 9; // add back dots on rearranged path (for poster)

     // Define input and output files.

     if ( !Member(vec<int>{1,2,3,7,8},display) )
         FatalErr("use DISPLAY=1,2,3,7,8 for this example");

     String adir, odir, tdir;
     SetTenxDirs( N, adir, odir, tdir );

     String dir  = DIR;         // output directories for SVG
     String obase = "xxx-" + ToString(N);  // use N to distinguish for caching
     String fout = "out.svg";

     // Get info.

     Bool cached = CACHED;
     if ( !cached )
     {    SystemSucceed( "NhoodInfo DIR_IN=" + adir +
               " OP=" + DIR +
               " O=" + obase + " S=" + SEEDS + " D=1 "
               "SHOW_ALIGN=False SHOW_CN=False PURPLE1=True ASPECT=0.2 SCALE=3 "
               "SVG=True" );
          // Take the pdf and change all colors to black, delete penwidth entry.
          // Reduce node width/height from 0.3 to 0.2.
          // Lower arrowsize from 3 to 2.
          // dot -Tsvg -o xxx.svg xxx.dot
          SystemSucceed( "EdgeInfo10X N=" + ToString(N) + " E=" + EDGES +
               " MODE=2 POS_REL=True NH=True > " + dir + "/" + obase );    }

     // Control.

     Bool print_text = False;

     // Output.

     Ofstream( sout, dir + "/" + fout );

     // Parse svg file.

     vec<String> head, tail;
     vec<svg_group> groups;
     ParseSvgFile( dir + "/" + obase + ".svg", head, groups, tail, print_text );
     vec<svg_group_plus> groupsp( groups.size( ) );
     for ( int i = 0; i < groups.isize( ); i++ )
          groupsp[i] = groups[i];

     // Print head and groups.

     for ( int i = 0; i < head.isize( ); i++ )
          sout << head[i] << endl;

     // Get points.

     int passes = 2;
     vec<vec<int>> point_edges {
//         {10258506,10269394},
//         {10254768,10265790}
         {10258506,10265790},
         {10254768,10269394}
     };
     vec<vec<int>> path_edges {
         {5246658,10258506,5332247,10265790,5412188},
         {5246658,10254768,5332247,10269394,5412188}
     };

     vec< vec< triple<int,int,double> > > points(4);
     vec< triple<int,int,double> > all_points;
     fast_ifstream pin( dir + "/" + obase );
     String line;
     while(1)
     {    getline( pin, line );
          if ( pin.fail( ) ) break;
          if ( !line.Contains( ":" ) ) continue;
          line.GlobalReplaceBy( ":", " " );
          istringstream iline( line.c_str( ) );
          int bc, e;
          double pos;
          iline >> bc >> e >> pos;
          all_points.push( bc, e, pos );    }
     for ( int pass = 1; pass <= passes; pass++ )
     {
          cout << "first PASS=" << pass << endl;
          for ( int r = 0; r < all_points.isize( ); r++ )
          {    int s;
               for ( s = r + 1; s < all_points.isize( ); s++ )
                    if ( all_points[s].first != all_points[r].first ) break;
               vec<int> es;
               for ( int t = r; t < s; t++ )
                    es.push_back( all_points[t].second );
               UniqueSort(es);
//               cout << "point_edges: " << printSeq(point_edges[pass-1]) << endl;
//               cout << "es: " << printSeq(es) << endl;
               Bool OK = True;
               if (!Meet(es, point_edges[pass-1])) OK=False;
#if 0
               if ( pass == 1 )
               {    if ( !Meet( es, {371196} ) ) OK = False;
                    if ( !Meet( es, {882135,7873161} ) ) OK = False;    }
               if ( pass == 2 )
               {    if ( !Meet( es, {8486660} ) ) OK = False;    }
               if ( pass == 3 )
               {    if ( !Meet( es, {7840986,8520273} ) ) OK = False;    }
               if ( pass == 4 )
               {    if ( !Meet( es, {371196} ) ) OK = False;
                    if ( !Meet( es, {371197,418206,418207,6261811} ) ) 
                         OK = False;    }
#endif
               if (OK)
               {    for ( int t = r; t < s; t++ )
                         points[pass-1].push_back( all_points[t] );    }
               r = s - 1;
          }
//          cout << "points: " << printSeq(points[pass-1]) << endl;
     }
               
     // Demonstrate adding a superpath.

     for ( int mp = 0; mp < passes; mp++ )
     {
          cout << "second PASS=" << mp+1 << endl;
          if ( display == 1 ) continue;
          if ( display == 2 && mp != 0 ) continue;
          if ( display == 3 && mp != 0 ) continue;
          if ( display == 4 && mp != 0 && mp != 2 ) continue;
          if ( display == 5 && mp != 0 && mp != 1 ) continue;
          if ( display == 6 && mp != 0 && mp != 3 ) continue;
          vec<int> es;
          es = path_edges[mp];
#if 0
          if ( mp == 0 ) es = {371196,8471468,418208,6249730,7873161,882135};
          if ( mp == 1 ) es = {576791,8486660,418208,6249730,7873161,882135};
          if ( mp == 2 ) es = {576791,7840986,6249730,8520273,882135};
          if ( mp == 3 ) es = {371196,371197,418207,6261811};
#endif

          pcb_path path;
          for ( int i = 0; i < es.isize( ); i++ )
          {    int e = es[i]; 
               cout << i << ": " << e << endl;
               if ( i == 0 ) path = EdgePath( groupsp, e );
               else path = Cat( path, EdgePath( groupsp, e ) );    }

          double delta = 0;
          if ( mp == 0 ) delta = 15;
          if ( mp == 1 ) delta = -18;
          if ( mp == 2 ) delta = 36;
          if ( mp == 3 ) delta = -18;

          const int dot_radius = 6;
          const double line_width = 18.0;
          double line_opacity = 1.0;
          String color;

          vec<vec<int>> C(4);
          C[0] = {100,85,73}; // peachpuff
          C[1] = {83,83,83};  // lightgray
          C[2] = {68,85,90};  // lightblue
          C[3] = {90,90,98};  // lavender
          color = SvgColor( C[mp] );

          vec<vec<int>> colors;

          const double min_dist = 8.0;
          const double min_dist_path = 25.0;
          const int ncolors = 100;
          const int min_color = 120;

          GetPalette( ncolors, colors, C, min_dist_path, min_dist, min_color );

          for ( int j = 0; j < path.isize( ); j++ )
               path[j].second += delta;

          svg_cpath x;
          x.SetColor( color );
          x.SetStrokeWidth( line_width );
          x.SetStrokeOpacity( line_opacity );
          x.SetCoords( path );
          sout << x.ToString( );

          if ( display == 8 ) continue;
          if ( display == 9 && mp > 0 ) continue;

          vec<int> bcs;
          for ( int pi = 0; pi < points[mp].isize( ); pi++ )
               bcs.push_back( points[mp][pi].first );
          UniqueSort(bcs);
          auto ccopy = colors;
          while ( bcs.size() > colors.size() ) {
              colors.append(ccopy);
          }
#if 0
          if ( bcs.size( ) > colors.size( ) ) 
          {    cout << "Need " << bcs.size( ) - colors.size( ) 
                    << " more colors." << endl;
               Scram(1);    }
#endif
          vec<String> X;
          for ( int pi = 0; pi < points[mp].isize( ); pi++ )
          {    int bc = points[mp][pi].first, e = points[mp][pi].second;
               int bp = BinPosition( bcs, bc );
               if ( display == 2 && bp != 8 ) continue;
               if ( display == 2 ) PRINT(bc);
               double d = points[mp][pi].third, x, y;
               PointLoc( groupsp, e, d, x, y );
               svg_ellipse el;
               el.SetColor( SvgColor( colors[bp] ) );
               el.SetRadius( dot_radius );
               el.SetPos( x, y + delta );
               X.push_back( el.ToString( ) );    }
          vec<int> shuffle;
          Shuffle( X.size( ), shuffle );
          for ( int i = 0; i < X.isize( ); i++ )
               sout << X[ shuffle[i] ];     }

     // Output tail.

     for ( int i = 0; i < groups.isize( ); i++ )
          sout << groups[i].ToString( );

     for ( int i = 0; i < tail.isize( ); i++ )
          sout << "[" << i+1 << "] " << tail[i] << endl;    }
