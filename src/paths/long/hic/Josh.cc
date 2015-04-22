///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Annotate map file to show LACHESIS or OLGA ordering.  Hardcoded.

#include "FastIfstream.h"
#include "MainTools.h"
#include "TokenizeString.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/Lines.h"

int main( )
{
     RunTime( );

     String target = "OLGA";
     // String target = "JOSH";

     // Load assembly.

     String dir = "/wga/scr4/jaffe/GapToy/49962/a.fin2";
     HyperBasevectorX hb;
     BinaryReader::readFile( dir + "/a.s.hbx", &hb );
     vec<int> inv;
     BinaryReader::readFile( dir + "/a.s.inv", &inv );
     vec<vec<vec<vec<int>>>> lines;
     BinaryReader::readFile( dir + "/a.s.lines", &lines );
     vec<int> tol;
     GetTol( hb, lines, tol );
     vec<int> INV( lines.size( ) );
     for ( int i = 0; i < lines.isize( ); i++ )
          INV[i] = tol[ inv[ lines[i][0][0][0] ] ];

     // Load Josh O/O, canonicalizing order.

     vec< pair< int, vec<int> > > J;
     String line;
     if ( target == "JOSH" )
     {    String josh = "/wga/scr4/tsharpe/lachesis/lachesisout/main_results";
          vec<String> all = AllFiles(josh);
          for ( int a = 0; a < all.isize( ); a++ )
          {    if ( !all[a].Contains( "group", 0 ) ) continue;
               int g = all[a].Between( "group", "." ).Int( );
               fast_ifstream in( josh + "/" + all[a] );
               vec<String> x;
               vec<int> tig;
               while(1)
               {    getline( in, line );
                    if ( in.fail( ) ) break;
                    if ( line.Contains( "#", 0 ) ) continue;
                    Tokenize( line, '\t', x );
                    int l = x[1].After( ":" ).Int( );
                    Bool rc = ( x[2] == "1" );
                    if (rc) l = INV[l];
                    tig.push_back(l);    }
               J.push( g, tig );    }    }
     
     // Load Olga O/O, canonicalizing order.

     if ( target == "OLGA" )
     {    fast_ifstream trans( dir + "/contig_id_translation" );
          vec<String> x;
          vec<int> to_line( 100000, -1 );
          while(1)
          {    getline( trans, line );
               if ( trans.fail( ) ) break;
               if ( line.Contains( "#", 0 ) ) continue;
               Tokenize( line, ' ', x );
               to_line[ x[1].Int( ) ] = x[0].Int( );    }
          String olga = dir + "/final.scaffold.list.txt";
          fast_ifstream in(olga);
          int count = 0;
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               Tokenize( line, ' ', x );
               vec<int> tig;
               for ( int i = 0; i < x.isize( ); i++ )
               {    int l = to_line[ Abs( x[i].Int( ) ) ];
                    if ( x[i].Int( ) < 0 ) l = INV[l];
                    tig.push_back(l);    }
               J.push( count, tig );
               count++;    }    }

     // Build index to J.

     vec< triple<int,int,Bool> > jind( lines.size( ), make_triple(-1,-1,False) );
     for ( int i = 0; i < J.isize( ); i++ )
     for ( int m = 0; m < J[i].second.isize( ); m++ )
     {    jind[ J[i].second[m] ] = make_triple( J[i].first, m, False );
          jind[ INV[J[i].second[m]] ] = make_triple( J[i].first, m, True );    }

     // Go through map.

     vec<String> lmap;
     vec<int> lmap_line;
     fast_ifstream zin( dir + "/a.s.lines.map" );
     int ll = -1;
     while(1)
     {    getline( zin, line );
          if ( zin.fail( ) ) break;
          String label;
          Bool flag = False;
          if ( !line.Contains( "line" ) ) 
          {    label = "";
               ll = -1;    }
          else
          {    int l = line.Between( "line[", "]" ).Int( );
               if ( jind[l].first >= 0 )
               {    label = ( jind[l].third ? "-" : "+" );
                    label += ToString( jind[l].first );
                    label += "." + ToString( jind[l].second );    }
               Bool ok = False;
               if ( ll < 0 ) ok = True;
               else if ( label == "" ) ok = True;
               else if ( jind[l].first == jind[ll].first 
                    && jind[l].second == jind[ll].second + 1
                    && !jind[l].third && !jind[ll].third )
               {    ok = True;    }
               else if ( jind[l].first == jind[ll].first 
                    && jind[l].second == jind[ll].second - 1
                    && jind[l].third && jind[ll].third )
               {    ok = True;    }
               if ( !ok ) flag = True;
               if ( label != "" ) ll = l;    }
          while ( label.isize( ) < 10 ) label += " ";
          while ( line.size( ) < 80 ) line += " ";
          cout << label << line;
          if (flag) cout << "  **********";
          cout << endl;    }    }
