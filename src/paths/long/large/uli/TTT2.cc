
#include "MainTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/Lines.h"

int main( )
{    RunTime( );
     String dir = "/wga/scr4/jaffe/GapToy/40";
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

     vec<int64_t> subsam_starts;
     BinaryReader::readFile( dir + "/subsam.starts", &subsam_starts );

     ReadPathVec paths( dir + "/a.final/a.paths" );

     subsam_starts.push_back( paths.size( ) );

     vec< vec<int> > delta( 9, vec<int>(9,0) );

     cout << "\ndifferentials:" << endl;
     for ( int e = 0; e < hb.E( ); e++ )
     {    int n = hb.Bases(e);
          int v = hb.ToLeft(e), w = hb.ToRight(e);
          if ( inv[e] < e ) continue;

          const int minc = 10;

          Bool print = False;
          vec<int64_t> ids;
          for ( int j1 = 0; j1 < 3; j1++ )
          for ( int j2 = 0; j2 < 3; j2++ )
          {    if ( count[j1][e] >= 2 && count[j2][e] >= 10 * count[j1][e] )
               {    for ( int64_t i = subsam_starts[j1]; 
                         i < subsam_starts[j1+1]; i++ )
                    for ( int l = 0; (int) l < (int) paths[i].size( ); l++ )
                    {    if ( paths[i][l] == e || paths[i][l] == inv[e] ) 
                              ids.push_back(i);    }
                    print = True;    }    }
          if (print)
          {    UniqueSort(ids);
               cout << "e = " << e << ", count = "
                    << count[0][e] << "," << count[1][e] << ","
                    << count[2][e] << ", ids = " << printSeq(ids) << endl;    }

          continue; // *************************************************************

          if ( count[0][e] > 0 ) continue;
          if ( count[1][e] < minc ) continue;
          if ( count[2][e] < minc ) continue;
          if ( count[3][e] > 0 ) continue;
          if ( count[4][e] < minc ) continue;
          if ( count[5][e] > 0 ) continue;
          if ( count[6][e] < minc ) continue;
          if ( count[7][e] > 0 ) continue;

          // Exclude edges starting or ending in long homopolymers.

          if ( hb.O(e).ToString( ).Contains( AN, 0 ) ) continue;
          if ( hb.O(e).ToString( ).Contains( AN, -1 ) ) continue;
          if ( hb.O(e).ToString( ).Contains( TN, 0 ) ) continue;
          if ( hb.O(e).ToString( ).Contains( TN, -1 ) ) continue;

          // Exclude dead ends.

          if ( hb.From(w).size( ) == 0 ) continue;
          if ( hb.To(v).size( ) == 0 ) continue;

          cout << "[" << ++C << "] ";
          PRINT(e);    }

     Scram(0);    }
