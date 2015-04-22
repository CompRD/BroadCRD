///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// FindCrossovers.  Find putative cancer rearrangements.  Hardwired now for
// certain samples.

#include "Intvector.h"
#include "MainTools.h"
#include "Qualvector.h"
#include "TokenizeString.h"
#include "VecUtilities.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(SAMPLE, 
          "NA12878 or cancer1 or normal1 or cancer2 or normal2");
     EndCommandArguments;

     String work_dir, FIN;
     if ( SAMPLE == "NA12878" ) 
     {    // work_dir = "/wga/scr4/jaffe/GapToy/49294b";
          work_dir = "/wga/scr4/jaffe/GapToy/49475.2";
          FIN = "fin";    }
     else if ( SAMPLE == "cancer1" ) 
     {    work_dir = "/wga/scr4/jaffe/GapToy/49264.HCC1954";
          FIN = "fin2";    }
     else if ( SAMPLE == "normal1" )
     {    work_dir = "/wga/scr4/jaffe/GapToy/49308.HCC1954BL";
          FIN = "fin2";    }
     else if ( SAMPLE == "cancer2" )
     {    work_dir = "/wga/scr4/jaffe/GapToy/49441.HCC1143";
          FIN = "fin";    }
     else if ( SAMPLE == "normal2" )
     {    work_dir = "/wga/scr4/jaffe/GapToy/49355.HCC1143BL";
          FIN = "fin";    }
     else
     {    cout << "Unknown sample.";
          Scram(1);    }
     String fin_dir = work_dir + "/a." + FIN;

     HyperBasevector hb;
     BinaryReader::readFile( fin_dir + "/a.hbv", &hb );
     vec<int> to_right;
     hb.ToRight(to_right);
     vec<int> inv;
     BinaryReader::readFile( fin_dir + "/a.inv", &inv );

     cout << ToStringAddCommas( hb.EdgeObjectCount( ) ) << " edges" << endl;

     vec< vec< pair<int,int> > > hits;
     BinaryReader::readFile( fin_dir + "/a.aligns", &hits );

     vec< vec<String> > places( hb.EdgeObjectCount( ) );
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
     {    int re = inv[e];
          int nhits = hits[e].size( );
          if ( re >= 0 ) nhits += hits[re].size( );
          if ( nhits == 0 ) continue;
          for ( int j = 0; j < hits[e].isize( ); j++ )
          {    places[e].push_back( "+" + ToString(hits[e][j].first) + "." 
                    + ToString(hits[e][j].second) );    }
          if ( re >= 0 )
          {    for ( int j = 0; j < hits[re].isize( ); j++ )
               {    places[e].push_back( "-" + ToString(hits[re][j].first)
                         + "." + ToString(hits[re][j].second) );    }    }    }


     //               e1'                                   e3
     // ---------------------------------> * --------------------------------->
     //                                    ^
     //                                   /  e2
     //               e1                 /                  e3'
     // ------------------------------> * ------------------------------------>

     int count = 0;
     for ( int e1 = 0; e1 < hb.EdgeObjectCount( ); e1++ )
     {    if ( !places[e1].solo( ) ) continue;
          int v = to_right[e1];
          if ( hb.From(v).size( ) != 2 ) continue;

          // Require that walking forward from v, the two edges do not arrive
          // at the same vertex within 4 kb, but do both extend at least 4 kb.

          const int min_apart = 100000;
          const int walk_dist = 4000;
          const int max_depth = 500; // NEED TO INCREASE!
          vec< vec< pair<int,int> > > down(2);
          for ( int i = 0; i < 2; i++ )
          {    int e = hb.EdgeObjectIndexByIndexFrom( v, i );
               down[i].push( hb.From(v)[i], hb.EdgeLengthKmers(e) );
               for ( int j = 0; j < down[i].isize( ); j++ )
               {    if ( down[i].isize( ) >= max_depth ) break;
                    if ( down[i][j].second >= walk_dist ) continue;
                    int z = down[i][j].first;
                    for ( int l = 0; l < hb.From(z).isize( ); l++ )
                    {    int f = hb.EdgeObjectIndexByIndexFrom( z, l );
                         down[i].push( hb.From(z)[l], down[i][j].second 
                              + hb.EdgeLengthKmers(f) );    }    }    }
          vec<vec<int>> dn(2);
          for ( int i = 0; i < 2; i++ )
          {    for ( int j = 0; j < down[i].isize( ); j++ )
                    dn[i].push_back( down[i][j].first );
               UniqueSort( dn[i] );    }
          if ( Meet( dn[0], dn[1] ) ) continue;
          int maxd = 0;
          Bool bad = False;
          for ( int i = 0; i < 2; i++ )
          {    int M = 0;
               for ( int j = 0; j < down[i].isize( ); j++ )
                    M = Max( maxd, down[i][j].second );
               if ( M < walk_dist ) bad = True;    }
          if (bad) continue;

          // Test as e2 each of the two edges starting from v

          for ( int i = 0; i < hb.From(v).isize( ); i++ )
          {    int e2 = hb.EdgeObjectIndexByIndexFrom( v, i );
               int e3p = hb.EdgeObjectIndexByIndexFrom( v, 1-i );
               if ( inv[e2] >= 0 && inv[e2] < e2 ) continue;
               int w = to_right[e2];
               if ( !hb.From(w).solo( ) ) continue;
               if ( hb.To(w).size( ) != 2 ) continue;
               int e1p = -1;
               for ( int j = 0; j < 2; j++ )
               {    if ( hb.EdgeObjectIndexByIndexTo( w, j ) != e2 )
                         e1p = hb.EdgeObjectIndexByIndexTo( w, j );    }

               // Require that walking backward from w, the two edges do not arrive
               // at the same vertex within 4 kb.

               vec< vec< pair<int,int> > > up(2);
               for ( int i = 0; i < 2; i++ )
               {    int e = hb.EdgeObjectIndexByIndexTo( w, i );
                    up[i].push( hb.To(w)[i], hb.EdgeLengthKmers(e) );
                    for ( int j = 0; j < up[i].isize( ); j++ )
                    {    if ( up[i].isize( ) >= max_depth ) break;
                         if ( up[i][j].second >= walk_dist ) continue;
                         int z = up[i][j].first;
                         for ( int l = 0; l < hb.To(z).isize( ); l++ )
                         {    int f = hb.EdgeObjectIndexByIndexTo( z, l );
                              up[i].push( hb.To(z)[l], up[i][j].second 
                                   + hb.EdgeLengthKmers(f) );    }    }    }
               vec<vec<int>> un(2);
               for ( int i = 0; i < 2; i++ )
               {    for ( int j = 0; j < up[i].isize( ); j++ )
                         un[i].push_back( up[i][j].first );
                    UniqueSort( un[i] );    }
               if ( Meet( un[0], un[1] ) ) continue;
               int maxu = 0;
               Bool bad = False;
               for ( int i = 0; i < 2; i++ )
               {    int M = 0;
                    for ( int j = 0; j < up[i].isize( ); j++ )
                         M = Max( maxd, up[i][j].second );
                    if ( M < walk_dist ) bad = True;    }
               if (bad) continue;

               // Go through one choice (stupid).

               for ( int j = 0; j < hb.From(w).isize( ); j++ )
               {    int e3 = hb.EdgeObjectIndexByIndexFrom( w, j );

                    /*
                    String top, bot;
                    if ( places[e3].solo( ) ) top = places[e3][0];
                    else if ( places[e1p].solo( ) ) top = places[e1p][0];
                    else continue;
                    if ( places[e1].solo( ) ) bot = places[e1][0];
                    else if ( places[e3p].solo( ) ) bot = places[e3p][0];
                    else continue;
                    if ( top.Before( "." ) == bot.Before( "." )
                         && Abs( top.After( "." ).Int( )
                               - bot.After( "." ).Int( ) ) < min_apart )
                    {    continue;    }
                    */

                    if ( !places[e3].solo( ) ) continue;
                    if ( places[e1][0].Before( "." ).Int( ) 
                         == places[e3][0].Before( "." ).Int( )
                         && Abs( places[e1][0].After( "." ).Int( )
                               - places[e3][0].After( "." ).Int( ) ) < min_apart )
                    {    continue;    }

                    // Get pids.
               
                    vec<int> ex = {e2};
                    if ( inv[e2] >= 0 ) ex.push_back( inv[e2] );
                    UniqueSort(ex);
                    VecULongVec x;
                    x.Read( fin_dir + "/a.paths.inv", ex );
                    vec<int> pids;
                    for ( int i = 0; i < (int) x.size( ); i++ )
                    for ( int j = 0; j < (int) x[i].size( ); j++ )
                         pids.push_back( x[i][j]/2 );
                    UniqueSort(pids);
                    int npids = pids.size( );

                    // Fetch paths.

                    vec<int> ids;
                    for ( int i = 0; i < npids; i++ )
                         ids.push_back( 2*pids[i], 2*pids[i] + 1 );
                    ReadPathVec p;
                    p.Read( fin_dir + "/a.paths", ids );

                    // Require both ends placed.

                    vec<Bool> to_delete( ids.size( ), False );
                    for ( int i = 0; i < npids; i++ )
                    {    if ( p[2*i].size( ) == 0 || p[2*i+1].size( ) == 0 )
                              to_delete[2*i] = to_delete[2*i+1] = True;    }
                    p.EraseIf(to_delete);

                    // Compute left and right support of x2.

                    vec<int> sleft, sright;
                    for ( int i = 0; i < (int) p.size( ); i++ )
                    {    vec<int> x;
                         for ( int j = 0; j < (int) p[i].size( ); j++ )
                              x.push_back( p[i][j] );
                         for ( int pass = 1; pass <= 2; pass++ )
                         {    if ( pass == 2 )
                              {    x.ReverseMe( );
                                   for ( int j = 0; j < x.isize( ); j++ )
                                        x[j] = inv[ x[j] ];    }
                              vec<int> l = { e1, e2 }, r = { e2, e3 };
                              if ( x.Contains(l) ) sleft.push_back(i/2);
                              if ( x.Contains(r) ) sright.push_back(i/2);    }    }
                    UniqueSort(sleft), UniqueSort(sright);
                    if ( sleft.isize( ) < 5 || sright.isize( ) < 5 ) continue;

                    cout << "\n[" << ++count << "] e1 = " << e1 << "[" 
                         << printSeq(places[e1]) << "], " << "e2 = " << e2 << "[" 
                         << printSeq(places[e2]) << "], " << "e3 = " << e3 << "[" 
                         << printSeq(places[e3]) << "]" << endl;    
                    PRINT2( sleft.size( ), sright.size( ) );    }    }    }
     Scram(0);    }
