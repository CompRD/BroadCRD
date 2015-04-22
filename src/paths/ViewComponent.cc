/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// ViewComponent: generate an html-viewable form of a single component from 
// a HyperKmerPath, along with ancillary truth data if available.
//
// The rendering of the dot file is crappy.  I don't know how to improve it.

#include "Basevector.h"
#include "MainTools.h"
#include "lookup/LookAlign.h"
#include "paths/AlignHyperKmerPath.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"

String Quote( const String& s )
{    return "\"" + s + "\"";    }

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_Int(K);
     CommandArgument_Bool_OrDefault(USE_TRUTH, False);
     CommandArgument_String_OrDefault(TMPDIR, "/tmp");
     CommandArgument_String(SUBDIR);
     CommandArgument_Int(COMPONENT);
     CommandArgument_Bool_OrDefault(LABEL_EDGES, False);
     EndCommandArguments;

     // Set up directories.

     String data_dir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;

     // Get HyperKmerPath component.

     HyperKmerPath h( HyperKmerPath( sub_dir + "/hyper" ), COMPONENT );

     // Align merged HyperKmerPath edges to reference.

     vec<look_align> aligns;
     vec< vec<int> > aligns_index;
     vecbasevector genome;
     KmerBaseBroker* kbb = 0;
     if (USE_TRUTH)
     {    kbb = new KmerBaseBroker( run_dir, K );
          AlignHyperKmerPath( h, kbb, data_dir + "/genome", TMPDIR, aligns,
               aligns_index );
          genome.ReadAll( data_dir + "/genome.fastb" );    }

     // Set up passes.  There are two, as we generate two htmls which one can toggle
     // back and forth between.

     for ( int pass = 1; pass <= 2; pass++ )
     {
          // Start html output.

          String dotfile = sub_dir + "/" + ToString(COMPONENT) 
               + ( pass == 1 ? "" : ".L" ) + ".dot";
          String htmlfile = sub_dir + "/" + ToString(COMPONENT)
               + ( pass == 1 ? "" : ".L" ) + ".html";
          String htmlfile2 = sub_dir + "/" + ToString(COMPONENT)
               + ( pass == 1 ? ".L" : "" ) + ".html";
          String pngfile = sub_dir + "/" + ToString(COMPONENT) 
               + ( pass == 1 ? "" : ".L" ) + ".png";
          String scaledpngfile = sub_dir + "/" + ToString(COMPONENT) 
               + ( pass == 1 ? "" : ".L" ) + ".scaled.png";
          Ofstream( html, htmlfile );
          html << "<!DOCTYPE HTML PUBLIC "
               << Quote( "-//IETF//DTD HTML//EN" ) << ">" << "\n";
          html << "<html>" << "\n";
          html << "<head>" << "\n";
          html << "<title>Component " + ToString(COMPONENT) + "</title>" << "\n";
          html << "</head>" << "\n";
          html << "<body bgcolor=" << Quote( "white" ) 
               << " text=" << Quote( "black" ) << ">" << "\n";
          html << "<font size=4>" << "\n";
          html << "<H2>Component " << ToString(COMPONENT) << "</H2>" << "\n";
          html << "<H5>PRE=" << PRE << "</H5>" << "\n";
          html << "<H5>DATA=" << DATA << "</H5>" << "\n";
          html << "<H5>RUN=" << RUN << "</H5>" << "\n";
          html << "<H5>SUBDIR=" << SUBDIR << "</H5>" << "\n";
          html << "<P><A HREF=" << Quote(htmlfile2) << "><IMG SRC=" 
               << Quote(scaledpngfile) << "</A></P>" << "\n";
          html << "<PRE>" << "\n";

          // Print text version of HyperKmerPath.

          if (USE_TRUTH)
          {    PrintAlignedHyperKmerPath( 
                    html, h, kbb, genome, aligns, aligns_index, False );    }
          else h.PrintSummaryPlus( html, 0, 0, 0, 0, 0, False );

          // Write dot form of merged graph.

          {    Ofstream( dot, dotfile );
               h.PrintSummaryDOT0w( dot, False, pass == 2, LABEL_EDGES );    }

          // Finish html output.  There is a horrible kludge here, because of shared
          // library differences.

          html << "</PRE>" << "\n";
          SystemSucceed( "cat " + dotfile + " | dot -Tpng > " + pngfile );
          int status = Csh( 
               "setenv LD_LIBRARY_PATH .:/util/lib:/usr/local/lib:/usr/lib:/lib; "
               "cat " + pngfile + " | pngtopnm | pnmscale 0.7 | pnmtopng "
               + "> " + scaledpngfile );
          if ( status != 0 )
          {    cout << "png scaling failed" << endl;
               exit(1);    }
          html << "</body>" << "\n";
          html << "</html>" << "\n";    }

     // Get absolute path of htmlfile.  This is bad, but I've forgotten the
     // right way to do it.

     String htmlfile = sub_dir + "/" + ToString(COMPONENT) + ".html";
     char htmlpath0[10000];
     realpath( htmlfile.c_str( ), htmlpath0 );
     String htmlpath(htmlpath0);
     if ( !htmlpath.Contains( "/wga", 0 ) )
     {    cout << "Warning: your path isn't in /wga, so you won't be able "
               << "to view\n" << "it using the Broad URL.\n\n";    }

     // Tell user where to go.  The Broad URL is specific to the WGA group at
     // the Broad Institute.  The path has to reside in /wga or it won't work.

     String dotfile = sub_dir + "/" + ToString(COMPONENT) + ".dot";
     cout << "Two files created (html and png).  If the html file is viewed with "
          << "Netscape or\nFirefox, click to enlarge it.  This is not necessary "
          << "with Safari.\n";
     cout << "-------------------------------------------------------------"
          << "-------------------\n";
     cout << "html file: " << htmlfile << "\n";    
     cout << "png file:  " << dotfile + ".scaled.png" << "\n";
     cout << "(There are also versions with \".L\" in them.)\n";
     cout << "Broad URL: " << "iwww.broad.mit.edu" << htmlpath << endl;
     cout << "(login info at /wga/dev1/.htpasswd/passwd-plaintext)\n\n";
     cout << "Click on image to toggle vertex-number visibility.\n\n";    }
