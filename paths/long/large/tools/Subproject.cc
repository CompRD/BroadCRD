///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Extract fastb and qualb files for the read pairs assigned to particular edges.

#include "Basevector.h"
#include "Intvector.h"
#include "MainTools.h"
#include "Qualvector.h"
#include "VecUtilities.h"
#include "feudal/PQVec.h"
#include "paths/long/ReadPath.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault_Doc(DIR, ".", 
          "looks for DIR/a.{fastb,inv,paths,paths.inv}");
     CommandArgument_String_Doc(E, "comma-separated list of edge ids or @fn "
          "where fn is a file having one edge id per line");
     CommandArgument_String_Doc(OUT_HEAD, "dump reads to OUT_HEAD.{fastb,qualb}");
     EndCommandArguments;

     // Parse edges.

     if ( !E.Contains( "@" ) && !E.Contains( "{" ) ) E = "{" + E + "}";
     vec<int> e;
     ParseIntSet( E, e );
         
     // Load inversion.

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

     // Dump reads.

     vec<int> ids;
     for ( int i = 0; i < npids; i++ )
          ids.push_back( 2*pids[i], 2*pids[i] + 1 );
     vecbasevector bases;
     VecPQVec quals;
     {    bases.Read( DIR + "/../data/frag_reads_orig.fastb", ids );
          quals.Read( DIR + "/../data/frag_reads_orig.qualp", ids );
          bases.WriteAll( OUT_HEAD + ".fastb" );
          quals.WriteAll( OUT_HEAD + ".qualp" );    }    }
