
#include "MainTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/Lines.h"

int main( )
{    RunTime( );
     String dir = "/wga/scr4/jaffe/GapToy/60";
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

     cout << "\ndifferentials:" << endl;
     for ( int e = 0; e < hb.E( ); e++ )
     {    int n = hb.Bases(e);
          int v = hb.ToLeft(e), w = hb.ToRight(e);
          if ( inv[e] < e ) continue;

          const int minc = 10;

          /*
          if ( ( count[0][e] == 0 && count[1][e] >= minc )
               || ( count[0][e] >= minc && count[1][e] == 0 ) )
          {    PRINT(e);    }
          */

          if ( ( count[0][e] <= 1 && count[1][e] >= minc * Max( 1, count[0][e] ) )
               || ( count[0][e] >= minc * Max( 1, count[1][e] )
                    && count[1][e] <= 1 ) )
          cout << "e = " << e << ", counts = " << count[0][e]
               << ", " << count[1][e] << endl;

          continue;

          // Exclude edges starting or ending in long homopolymers.

          if ( hb.O(e).ToString( ).Contains( AN, 0 ) ) continue;
          if ( hb.O(e).ToString( ).Contains( AN, -1 ) ) continue;
          if ( hb.O(e).ToString( ).Contains( TN, 0 ) ) continue;
          if ( hb.O(e).ToString( ).Contains( TN, -1 ) ) continue;

          if ( count[0][e] > 0 ) continue;
          if ( count[1][e] < minc ) continue;
          if ( count[2][e] < minc ) continue;
          if ( count[3][e] > 0 ) continue;
          if ( count[4][e] < minc ) continue;
          if ( count[5][e] > 0 ) continue;
          if ( count[6][e] < minc ) continue;
          if ( count[7][e] > 0 ) continue;

          // Exclude dead ends.

          if ( hb.From(w).size( ) == 0 ) continue;
          if ( hb.To(v).size( ) == 0 ) continue;

          cout << "[" << ++C << "] ";
          PRINT(e);    }

     Scram(0);    }
