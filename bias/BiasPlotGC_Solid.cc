/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

const char *DOC =

"For each read perfect to N bases, associate the GC content "
"of the W-base window on the reference, starting at the beginning of the read.  "
"For each possible GC content, find the number of windows on the reference "
"having that GC content.  Then to each such GC content, associate "
"ratio (# of reads) / (# of windows on reference).  Normalize so that GC = 50% "
"has height 1."

;

#include "Basevector.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "graphics/BasicGraphics.h"
#include "lookup/PerfectCount.h"
#include "math/Functions.h"
#include "solexa/SolexaTools.h"
#include "solexa/SolexaMetrics.h"
#include <math.h>

#define BRACKFAIL                             \
{    cout << "Illegal bracketing in HEAD.\n"; \
     cout << "Abort.\n";                      \
     exit(1);    }

int main( int argc, char *argv[] ) 
{
     RunTime( );

     BeginCommandArguments;
     CommandDoc(DOC);
     CommandArgument_String_Doc(HEAD, 
          "flowcell.lane, or date_flowcell.lane, or list of such, separated by "
          "commas to keep separate or plus signs to merge; expressions like "
          "4321.{1,2,3} or 4321.{1+2} are formally expanded; at most 8 lanes in "
          "total are allowed; do not place entire list in brackets { }" );
     CommandArgument_String_Doc(GC_NAME, "Filename of GC vector.");
     CommandArgument_String_Doc(GC_REF_NAME, "Filename of GC_ref vector.");
     CommandArgument_Bool_OrDefault_Doc(ADDMETRIC, False,
          "If True, add gc window flatness to the metrics table." );
     CommandArgument_Bool_OrDefault(WRITE, True);
     CommandArgument_Bool_OrDefault_Doc(DOGRAPH, True,
          "If True, output graph." );
     CommandArgument_Bool_OrDefault_Doc(PRINTPOINTS, False,
          "If True, print to stdout the plot points." );
     CommandArgument_Int_OrDefault_Doc(W, 50, "window size");
     CommandArgument_Double_OrDefault_Doc(MAX_SLOP, 0.05, "spread of 95% "
          "confidence interval must be less than this to plot point" );
     CommandArgument_String_OrDefault_Doc(OUT, "", "output file, may end with .png");
     EndCommandArguments;

     // Check arguments.

     ForceAssertEq( W % 2, 0 );
     if (DOGRAPH) CheckValidGraphicsSuffix( OUT, True );

     // Parse HEAD into HEADS.

     vec< vec<String> > HEADS;
     String HEAD2, CURHEAD;
     for ( int i = 0; i < HEAD.isize( ); i++ )
     {    if ( HEAD[i] != '{' )
          {    HEAD2 += HEAD[i];
               continue;    }
          int j;
          for ( j = i + 1; j < HEAD.isize( ); j++ )
               if ( HEAD[j] == '}' ) break;
          if ( i < 5 || HEAD[i-1] != '.' || j == HEAD.isize( ) ) BRACKFAIL
          String brack = HEAD.substr( i, j - i + 1 );
          vec<int> lanes;
          Bool plus = False;
          if ( brack.Contains( "+" ) )
          {    if ( brack.Contains( "," ) ) BRACKFAIL
               plus = True;
               brack.GlobalReplaceBy( "+", "," );    }
          ParseIntSet( brack, lanes );
          if ( lanes.empty( ) ) BRACKFAIL
          int fcdotlen = 5;
          if ( HEAD2.size( ) > 5 && HEAD2[ HEAD2.isize( ) - 6 ] == '_' )
          {    fcdotlen = 12;
               ForceAssertGe( HEAD2.isize( ), fcdotlen );    }
          String fcdot = HEAD2.substr( HEAD2.isize( ) - fcdotlen, fcdotlen );
          HEAD2.resize( HEAD2.isize( ) - fcdotlen );
          for ( int k = 0; k < lanes.isize( ); k++ )
          {    if ( k > 0 ) HEAD2 += ( plus ? '+' : ',' );
               HEAD2 += fcdot + ToString( lanes[k] );    }
          i = j;    }
     vec<String> HEADPILE;
     for ( int i = 0; i < HEAD2.isize( ); i++ )
     {    if ( HEAD2[i] == ',' )
          {    HEADPILE.push_back(CURHEAD), HEADS.push_back(HEADPILE);
               HEADPILE.clear( ), CURHEAD.clear( );    }
          else if ( HEAD2[i] == '+' )
          {    HEADPILE.push_back(CURHEAD), CURHEAD.clear( );    }
          else CURHEAD += HEAD2[i];    }
     HEADPILE.push_back(CURHEAD), HEADS.push_back(HEADPILE);
     for ( int i = 0; i < HEADS.isize( ); i++ )
          Sort( HEADS[i] );

     // Create compressed version of each HEADS[i], for title.

     vec<String> TITLES;
     for ( int hi = 0; hi < HEADS.isize( ); hi++ )
     {    vec<String> H = HEADS[hi];
          String TITLE;
          for ( int i = 0; i < H.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < H.isize( ); j++ )
                    if ( H[j].Before( "." ) != H[i].Before( "." ) ) break;
               if ( i > 0 ) TITLE += "+";
               if ( j - i == 1 ) TITLE += H[i];
               else
               {    vec<int> lanes;
                    for ( int k = i; k < j; k++ )
                         lanes.push_back( H[k].After( "." ).Int( ) );
                    ForceAssert( lanes.UniqueOrdered( ) );
                    TITLE += H[i].Before( "." ) + ".{";
                    if ( lanes.isize( ) == lanes.back( ) - lanes.front( ) + 1
                         && lanes.size( ) > 3 )
                    {    TITLE += ToString( lanes.front( ) ) + "+...+"
                              + ToString( lanes.back( ) );    }
                    else
                    {    for ( int k = i; k < j; k++ )
                         {    if ( k > i ) TITLE += "+";
                              TITLE += ToString( lanes[k-i] );    }    }
                    TITLE += "}";    }
               i = j - 1;    }
          TITLES.push_back(TITLE);    }

     // Set up to plot points.  We define a fixed list of eight colors.

     vec<graphics_primitive> points;
     vec<color> colors;

     colors.push_back( red, blue, green );
     colors.push( 1, 0.25, 1 ),           colors.push( 1, 0.85, 0 );
     colors.push( 0, 0.8, 0.8 ),          colors.push( 0.65, 0.65, 0.65 );
     colors.push( 0.7, 0.3, 0.4 );

     // Go through the lane groups.

     ForceAssertLe( HEADS.size( ), colors.size( ) );
     for ( int hi = 0; hi < HEADS.isize( ); hi++ )
     {    
          // Get list of lanes in this group.

          vec<String> HEAD = HEADS[hi];

          // Create HEADXS = full paths to HEADS.

          vec<String> HEADXS;
          String HEADX = "{";
          for ( int i = 0; i < HEAD.isize( ); i++ )
          {    if ( i > 0 ) HEADX += ",";
               HEADX += HEAD[i];    }
          HEADX += "}";
          ExpandHead( HEADX, HEADXS );

          // Print out stdout header
          if (PRINTPOINTS) {
             cout << "Head #\t\tPercent\t\tGC Bias\n";   
          }



          // Load GC vectors
          // vec<int> GC( W+1, 0 ), GC_ref( W+1, 0 );
          vec<int> GC( W+1 ), GC_ref( W+1 );
          Ifstream(gc_in, GC_NAME);
          Ifstream(gc_ref_in, GC_REF_NAME);
          for ( int i = 0; i < W+1; i++ ) {
               gc_in >> GC[i];  
               gc_ref_in >> GC_ref[i];  
          }
          gc_in.close();
          gc_ref_in.close();

          // Generate GC representation ratios.
          vec<double> GC_ratio(W+1);
          for ( int i = 0; i <= W; i++ )
          {    if ( GC_ref[i] > 0 )
                        GC_ratio[i] = double( GC[i] ) / double( GC_ref[i] );    }
          if ( GC_ref[W/2] == 0 )
          {    cout << "Your genome is too weird.  "
                    << "It has no windows of GC content 50%.\n" << "Giving up.\n";
                    exit(1);    }
          if ( GC[W/2] == 0 )
          {    cout << "Your reads are too weird.  "
                    << "None found for window GC content 50%.\n" << "Giving up.\n";
          exit(1);    }
          double mid = GC_ratio[W/2];
              for ( int i = 0; i <= W; i++ )
                 if ( GC_ref[i] > 0 ) GC_ratio[i] /= mid;


          // Generate points to plot.

          double p = 0.95;
          double z = InverseNormalCDF( 1.0 - (1.0-p)/2.0 );
          vec<double> x, y;
          for ( int i = 0; i <= W; i++ )
          {    double u = GC[i];
               if ( GC_ref[i] > 0 && GC[i] > 0 )
               {    double low = u - z * sqrt(u), high = u + z * sqrt(u);
                    low = Max( low, 0.0 );
                    if ( (high-low) / ( GC_ref[i] * mid ) < MAX_SLOP )
                    // if ( low > 0 && high/low <= 1.05 )
                    {    x.push_back( 100.0 * double(i)/double(W) );
                         y.push_back( GC_ratio[i] );    
                    }    
               
               }    
          }
          double POINTSIZE = 1.0;
          points.push_back( SetColor( colors[hi] ) );

          if (x.isize() == 0 || x.isize() == 1) {
            cout << "No or just one perfect match found so the metric makes no sense.  Exiting.\n";
            exit(1);
          } else {
            double biasGCValL2 = 0.0;
            double h, fx1, fx0;
            double length = x[x.isize()-1] - x[0];
  
            for ( int i = 0; i < x.isize( ); i++ ) {
                 points.push_back( Point( x[i], y[i], POINTSIZE ) );    
                 if (i < x.isize()-1) {
                    h = x[i+1] - x[i];
                    fx0 = (1.0 - y[i])*(1.0 - y[i]);
                    fx1 = (1.0 - y[i+1])*(1.0 - y[i+1]);
                    biasGCValL2 += (h/2.0)*(fx0 + fx1);
                 }
                 if (PRINTPOINTS) {
                         cout << hi << "\t\t" << x[i] << "\t\t" << y[i] << "\n";
                 }
            }
            // L2 norm below
            biasGCValL2 = sqrt(biasGCValL2);
            // Normalize by area since the domain of the data varies.
            // Correct name of metric is "biasGCValL2Normalized".
            String bias_GC = ToString( 100.0 * biasGCValL2/length, 2 );
            PRINT(bias_GC);
            if (ADDMETRIC) {
              for ( int u = 0; u < HEADXS.isize( ); u++ ) {
                solexa_metric_db db(HEADXS[u] + ".metrics");
                db.SetValue("bias_GC", bias_GC);
                if (WRITE) db.WriteMetrics( HEADXS[u] + ".metrics" );
                else db.WriteMetrics(cout);
              }
            }
          }
     }




   if (DOGRAPH) {

     // Form title.

     const double TITLE_FONTSIZE = 15;
     vec<graphics_primitive> title;
     vec<String> heads0;
     for ( int i = 0; i < TITLES.isize( ); i++ )
          heads0.push_back( TITLES[i] + " " );
     colors.resize( HEADS.size( ) );
     title.push_back( RainbowTextCenter( heads0, colors, 100.0/2.0, 0, 0,
          TimesBold(TITLE_FONTSIZE) ) );

     // Generate graph.

     double LINEWIDTH = 0.5;
     double POSTSCALE = 0.2;
     points.push_back( SetColor(black), SetLineWidth(LINEWIDTH) );
     vec<graphics_primitive> x_axis = AxisX( 0, 100, 1.0, True, "%", 0.0 );
     x_axis.push_front( SetLineWidth(LINEWIDTH) );
     points.append( AxisY( 0, 100, 0, 2.3, True, "", 0.0 ) );
     vec<graphics_primitive> xlabel;
     xlabel.push_back( TextCenter( "GC content of " + ToString(W) + "bp window",
          100.0/2.0, 0, 0, TimesBold(11) ) );
     vec< vec<graphics_primitive> > stack;
     vec<double> heights;
     stack.push_back(xlabel),    heights.push_back(65);
     stack.push_back(x_axis),    heights.push_back(0);
     stack.push_back(points),    heights.push_back(200);
     stack.push_back(title);     heights.push_back(20);
     RenderGraphics( OUT, stack, heights, 1.0, 200, 50, True, POSTSCALE );    
     
  // End of if (DOGRAPH) {
  }     
     
}
