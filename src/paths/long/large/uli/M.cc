

// order  assembly  true   state   sample                    
// 0      mouse1    #1     0       Balb/c, male              
// 1      mouse4    #4b    1       son, responder            
// 2      mouse5    #5     1       son, responder            
// 3      mouse6    #6     0       son, nonresponder         
// 4      mousea    #4a    1       son, responder
// 5      mouseb    #7     0       son, nonresponder        
// 6      mousec    #3     1       ~B6 mother, responder
// 7      moused    #2     0       B6 father
// 8      mousee    DUP    1       reads reordered from mousea

//          0        1        2        3        4        5        6        7        8
// 0        0  2005649  2036713  2038518  2025782  2026224  2028392  2059196  2025786
// 1  2005649        0     4409     2117      122     4783    47907     4207      124
// 2  2036713     4409        0     3433     4441     3928    50235     3276     4441
// 3  2038518     2117     3433        0     2131     6538    50181     5653     2130
// 4  2025782      122     4441     2131        0     4747    52902     5336        0
// 5  2026224     4783     3928     6538     4747        0    52621     4826     4749
// 6  2028392    47907    50235    50181    52902    52621        0    38343    52889
// 7  2059196     4207     3276     5653     5336     4826    38343        0     5330
// 8  2025786      124     4441     2130        0     4749    52889     5330        0

#include "MainTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/Lines.h"

int main( )
{
     String dir = "/wga/scr4/jaffe/GapToy/mouse.1456abcde";
     HyperBasevectorX hb;
     BinaryReader::readFile( dir + "/a.final/a.hbx", &hb );
     vec<int> inv;
     BinaryReader::readFile( dir + "/a.final/a.inv", &inv );
     vec<vec<vec<vec<int>>>> lines;
     BinaryReader::readFile( dir + "/a.final/a.lines", &lines );
     vec<int> tol;
     GetTol( hb, lines, tol );
     vec<vec<int>> count;
     BinaryReader::readFile( dir + "/a.final/a.countsb", &count );
     String AN = String( 20, 'A' );
     String TN = String( 20, 'T' );
     int C = 0;

     vec< vec<int> > delta( 9, vec<int>(9,0) );

     int c1 = 0, c2 = 0;

     int c12 = 0, c13 = 0, c15 = 0, c23 = 0, c25 = 0, c35 = 0;

     for ( int e = 0; e < hb.E( ); e++ )
     {    int n = hb.Bases(e);
          int v = hb.ToLeft(e), w = hb.ToRight(e);
          if ( inv[e] < e ) continue;

          // Require bubbles.

          if ( hb.From(v).size( ) != 2 ) continue;
          if ( hb.From(v)[0] != w ) continue;
          if ( hb.From(v)[1] != w ) continue;
          if ( hb.To(w).size( ) != 2 ) continue;
          if ( v == w ) continue;

          int minc = 10;

          /*
          if ( count[0][e] > 0 ) continue;
          if ( count[1][e] < minc ) continue;
          if ( count[2][e] < minc ) continue;
          if ( count[3][e] > 0 ) continue;
          if ( count[4][e] < minc ) continue;
          if ( count[5][e] > 0 ) continue;
          if ( count[6][e] < minc ) continue;
          if ( count[7][e] > 0 ) continue;
          */

          // Exclude edges starting or ending in long homopolymers.

          if ( hb.O(e).ToString( ).Contains( AN, 0 ) ) continue;
          if ( hb.O(e).ToString( ).Contains( AN, -1 ) ) continue;
          if ( hb.O(e).ToString( ).Contains( TN, 0 ) ) continue;
          if ( hb.O(e).ToString( ).Contains( TN, -1 ) ) continue;

          // Exclude dead ends.

          if ( hb.From(w).size( ) == 0 ) continue;
          if ( hb.To(v).size( ) == 0 ) continue;


          /*
          if ( ( ( count[8][e] >= minc && count[4][e] == 0 )
               || ( count[4][e] >= minc && count[8][e] == 0 ) )
               && hb.Kmers(e) > 1000 )
          {    cout << "[" << ++C << "] ";
               PRINT(e);    }
          continue;
          */

          /*
          minc = 20;
          if ( hb.Kmers(e) < 200 ) continue;
          if ( count[0][e] >= minc )
          {    if ( count[1][e] >= minc && count[7][e] == 0 ) c1++;
               if ( count[7][e] >= minc && count[1][e] == 0 ) c2++;    }
          continue;
          */

          // Look for strong differences between some of the mice.

          /*
          minc = 20;
          if ( hb.Kmers(e) < 200 ) continue;
          for ( int i1 = 1; i1 < 8; i1++ )
          for ( int i2 = 1; i2 < 8; i2++ )
          {    if ( i1 == 4 || i2 == 4 ) continue;
               if ( i1 == 6 || i2 == 6 ) continue;
               if ( count[i1][e] >= minc && count[i2][e] == 0 )
               {    PRINT3( e, i1, i2 );    }    }
          continue;
          */

          int m1 = count[1][e], m2 = count[2][e], m3 = count[3][e], m5 = count[5][e];
          if ( m1 >= minc && m2 >= minc && m3 == 0 && m5 == 0 ) c12++;
          if ( m1 >= minc && m3 >= minc && m2 == 0 && m5 == 0 ) c13++;
          if ( m1 >= minc && m5 >= minc && m2 == 0 && m3 == 0 ) c15++;
          if ( m2 >= minc && m3 >= minc && m1 == 0 && m5 == 0 ) c23++;
          if ( m2 >= minc && m5 >= minc && m1 == 0 && m3 == 0 ) c25++;
          if ( m3 >= minc && m5 >= minc && m1 == 0 && m2 == 0 ) c35++;
          continue;

          // Make yet another table.

          /*
          if ( hb.Kmers(e) < 200 ) continue;
          for ( int i1 = 0; i1 < 9; i1++ )
          for ( int i2 = 0; i2 < 9; i2++ )
          if ( count[i1][e] >= minc && count[i2][e] == 0 )
          {    delta[i1][i2]++;    }
          continue;
          */

          // Make different table.

          /*
          if ( count[0][e] < minc ) continue; // conditional on Balb/c.
          if ( hb.Kmers(e) < 200 ) continue;
          for ( int i1 = 0; i1 < 9; i1++ )
          for ( int i2 = 0; i2 < 9; i2++ )
          if ( count[i1][e] >= minc && count[i2][e] == 0 )
          {    delta[i1][i2]++;    }
          continue;
          */

          // Make table.

          /*
          if ( hb.Kmers(e) < 200 ) continue;
          for ( int i1 = 0; i1 < 9; i1++ )
          for ( int i2 = 0; i2 < 9; i2++ )
          if ( ( count[i1][e] >= minc && count[i2][e] == 0 )
               || ( count[i2][e] >= minc && count[i1][e] == 0 ) )
          {    delta[i1][i2]++;    }
          continue;
          */


          /*
          if ( hb.From(v).size( ) != 2 ) continue;
          if ( hb.From(v)[0] != w ) continue;
          if ( hb.From(v)[1] != w ) continue;
          if ( hb.To(w).size( ) != 2 ) continue;
          if ( v == w ) continue;

          int ep = -1;
          if ( hb.IFrom(v,0) != e ) ep = hb.IFrom(v,0);
          if ( hb.IFrom(v,1) != e ) ep = hb.IFrom(v,1);
          if ( count[1][ep] > 0 ) continue;
          if ( count[2][ep] > 0 ) continue;
          if ( count[4][ep] > 0 ) continue;
          if ( count[6][ep] > 0 ) continue;
          */

          cout << "[" << ++C << "] ";
          PRINT(e);    }

     
     /*
     vec<vec<String>> rows;
     vec<String> row = {""};
     for ( int i = 0; i < 9; i++ )
          row.push_back( ToString(i) );
     rows.push_back(row);
     for ( int i1 = 0; i1 < 9; i1++ )
     {    vec<String> row;
          row.push_back( ToString(i1) );
          for ( int i2 = 0; i2 < 9; i2++ )
               row.push_back( ToString( delta[i1][i2] ) );
          rows.push_back(row);    }
     PrintTabular( cout, rows, 2, "rrrrrrrrrr" );
     */

     // PRINT2( c1, c2 );

     PRINT6( c12, c13, c15, c23, c25, c35 );

     Scram(0);    }
