///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// TestDegloop.  Run Degloop on particular edges.

#include "Basevector.h"
#include "Intvector.h"
#include "MainTools.h"
#include "VecUtilities.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault_Doc(DIR, ".", 
          "looks for DIR/a.{fastb,inv,paths,paths.inv}");
     CommandArgument_String_Doc(E, "comma-separated list of edges to use");
     CommandArgument_Int_OrDefault_Doc(VERBOSITY, 2, "or 3");
     CommandArgument_Int_Doc(MODE, "1 or 2");
     CommandArgument_Int_OrDefault(MIN_DIST, 3.0);
     EndCommandArguments;
         
     // Parse E.

     vec<int> e;
     ParseIntSet( E, e );
     if ( e.size( ) != 2 )
     {    cout << "Please supply two edges." << endl;
          Scram(1);    }

     // Load assembly.

     HyperBasevectorX hb;
     BinaryReader::readFile( DIR + "/a.hbx", &hb );
     if ( hb.ToLeft( e[0] ) != hb.ToLeft( e[1] )
          && hb.ToRight( e[0] ) != hb.ToRight( e[1] ) )
     {    cout << "Those two edges don't share a vertex." << endl;
          Scram(1);    }
     vec<int> inv;
     BinaryReader::readFile( DIR + "/a.inv", &inv );

     // Look up edges.

     for ( int i = 0; i < e.isize( ); i++ )
     {    if ( e[i] >= inv.isize( ) )
          {    cout << e[i] << " is not a valid edge id." << endl;
               Scram(1);    }    }
     int ne = e.size( );
     vec<int> e2(e);
     for ( int i = 0; i < ne; i++ )
          if ( inv[ e[i] ] >= 0 ) e2.push_back( inv[ e[i] ] );
     UniqueSort(e2);
     VecULongVec x0;
     x0.Read( DIR + "/a.paths.inv", e2 );
     vec<vec<int>> x( x0.size( ) );
     for ( int i = 0; i < x.isize( ); i++ )
     for ( int j = 0; j < (int) x0[i].size( ); j++ )
          x[i].push_back( x0[i][j] );

     // Define pids.

     vec<int> pids;
     for ( int j = 0; j < e2.isize( ); j++ )
     for ( int i = 0; i < (int) x[j].size( ); i++ )
          pids.push_back( x[j][i]/2 );
     UniqueSort(pids);
     int npids = pids.size( );

     // Fetch paths and reads.

     vec<int> ids;
     for ( int i = 0; i < npids; i++ )
          ids.push_back( 2*pids[i], 2*pids[i] + 1 );
     ReadPathVec paths;
     paths.Read( DIR + "/a.paths", ids );
     vecbasevector bases;
     bases.Read( DIR + "/../data/frag_reads_orig.fastb", ids );
     VecPQVec quals;
     quals.Read( DIR + "/../data/frag_reads_orig.qualp", ids );

     // Set up paths_index.

     VecULongVec paths_index( hb.E( ) );
     for ( int i = 0; i < e2.isize( ); i++ )
     for ( int j = 0; j < (int) x0[i].size( ); j++ )
          paths_index[ e2[i] ].push_back( BinPosition( ids, x0[i][j] ) );

     // Degloop.

     vec<int> EDELS;
     if ( hb.ToLeft( e[0] ) == hb.ToLeft( e[1] ) )
     {    DegloopCore( MODE, hb, inv, paths, bases, quals, paths_index, 
               hb.ToLeft(e[0]), 1, MIN_DIST, EDELS, VERBOSITY, &ids );    }
     if ( hb.ToRight( e[0] ) == hb.ToRight( e[1] ) )
     {    DegloopCore( MODE, hb, inv, paths, bases, quals, paths_index, 
               hb.ToRight(e[0]), 2, MIN_DIST, EDELS, VERBOSITY, &ids );    }    }
