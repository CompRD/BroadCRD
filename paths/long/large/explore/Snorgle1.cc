///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// This is an experimental version of part of Snorgle phase 1.  Needs to be kept 
// in sync with it.

#include "Basevector.h"
#include "MainTools.h"
#include "Qualvector.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_Int_OrDefault_Doc(TRACE_PID, -1,
          "for phase 2, analyze just this PID; also follow in phase 1");
     CommandArgument_Int_OrDefault_Doc(TRACE_RID, -1,
          "sets TRACE_PID = TRACE_RID/2");
     CommandArgument_String_OrDefault_Doc(INSTANCE, "51400.newchem",
          "s or 1 or 2 or 3 or 4 or 5 or 51400.newchem");
     EndCommandArguments;

     if ( TRACE_RID >= 0 ) TRACE_PID = TRACE_RID/2;
     ForceAssertGe( TRACE_PID, 0 );

     String work_dir = "/wga/scr4/jaffe/GapToy/" + INSTANCE;

     int64_t id1 = TRACE_PID * 2;
     int64_t id2 = id1 + 1;
     vec<int64_t> ids = {id1,id2};

     cout << Date( ) << ": loading" << endl;
     ReadPathVec paths;
     paths.Read( work_dir + "/a.200/a.paths", ids );
     vecbasevector bases;
     bases.Read( work_dir + "/data/frag_reads_orig.fastb", ids );
     VecPQVec pquals;
     pquals.Read( work_dir + "/data/frag_reads_orig.qualp", ids );
     HyperBasevectorX hb;
     BinaryReader::readFile( work_dir + "/a.200/a.hbx", &hb );
     vec<int> inv;
     BinaryReader::readFile( work_dir + "/a.200/a.inv", &inv );
     
     // Heuristics.

     const int adaptn = 20;
     
     // Get distance to end for each vertex.

     vec<int> D;
     BinaryReader::readFile( work_dir + "/dist_to_end", &D );

     // Kill reads having adapter.
     
     double clock = WallClockTime( );
     cout << Date( ) << ": killing reads having adapter" << endl;
     String adapt 
          = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTATCTCGTATGCCGTCTTCTGCTTG";
     String adapt20 = adapt;
     adapt20.resize(adaptn);

     // Identify incompletely placed pairs, having some quality.  Ignore reads 
     // containing adapter.

     const int min_qual = 1000;
     const int max_prox = 500;
     int64_t pid = TRACE_PID;
     {    int64_t id1 = 2*pid, id2 = 2*pid+1;
          const ReadPath &p1 = paths[0], &p2 = paths[1];
          cout << "p1: " << p1.getOffset( ) << " : " << printSeq(p1) << endl;
          cout << "p2: " << p2.getOffset( ) << " : " << printSeq(p2) << endl;

          // Test proximity.

          int d1 = -1, d2 = -1;
          if ( p1.size( ) > 0 )
          {    int start1 = p1.getOffset( );
               int stop1 = start1 + bases[0].isize( );
               for ( int j = 1; j < (int) p1.size( ); j++ )
                    stop1 -= hb.Kmers( p1[j] );
               int d1a = hb.Bases( p1.back( ) ) - stop1
                    + D[ hb.ToRight( p1.back( ) ) ];
               int d1b = start1 + D[ hb.ToRight( inv[ p1.front( ) ] ) ];
               d1 = Min( d1a, d1b );
               PRINT3( d1, d1a, d1b );    }
          if ( p2.size( ) > 0 )
          {    int start2 = p2.getOffset( );
               int stop2 = start2 + bases[1].isize( );
               for ( int j = 1; j < (int) p2.size( ); j++ )
                    stop2 -= hb.Kmers( p2[j] );
               int d2a = hb.Bases( p2.back( ) ) - stop2
                    + D[ hb.ToRight( p2.back( ) ) ];
               int d2b = start2 + D[ hb.ToRight( inv[ p2.front( ) ] ) ];
               d2 = Min( d2a, d2b );
               PRINT3( d2, d2a, d2b );    }
          if ( d1 > max_prox && d2 > max_prox )
          {    cout << "distance to end = " << d1 << ", too far" << endl;
               cout << "distance to end = " << d2 << ", too far" << endl;
               return 0;    }
          if ( d1 > max_prox && d2 == -1 )
          {    cout << "distance to end = " << d1 << ", too far" << endl;
               return 0;    }
          if ( d2 > max_prox && d1 == -1 )
          {    cout << "distance to end = " << d2 << ", too far" << endl;
               return 0;    }

          cout << "past first check" << endl;

          /*
          if ( p1.size( ) > 0 && p2.size( ) > 0 )
          {    int e1 = p1.back( ), e2 = p2.back( );
               // if ( hb.From( hb.ToRight(e1) ).nonempty( ) ) return 0;
               // if ( hb.From( hb.ToRight(e2) ).nonempty( ) ) return 0;
               if ( inv[e1] == e2 ) return 0;    }
          cout << "past second check" << endl;
          */

          int qsum1 = 0, qsum2 = 0;
          qvec qv1, qv2;
          pquals[0].unpack( &qv1 ), pquals[1].unpack( &qv2 );
          for ( int j = 0; j < (int) qv1.size( ); j++ )
               if ( qv1[j] > 2 ) qsum1 += qv1[j];
          for ( int j = 0; j < (int) qv2.size( ); j++ )
               if ( qv2[j] > 2 ) qsum2 += qv2[j];

          if ( p1.size( ) == 0 && qsum1 < min_qual ) 
          {    cout << "read 1 unplaced and has qsum " << qsum1 << endl;
               return 0;    }
          if ( p2.size( ) == 0 && qsum2 < min_qual ) 
          {    cout << "read 2 unplaced and has qsum " << qsum2 << endl;
               return 0;    }
          cout << "past third check" << endl;

          Bool have_ad = False;
          String s1 = bases[0].ToString( );
          for ( int j = 0; j <= bases[0].isize( ) - adaptn; j++ )
          {    if ( s1.Contains( adapt20, j ) )
               {    have_ad = True;
                    break;    }    }
          if (have_ad) return 0;
          String s2 = bases[1].ToString( );
          for ( int j = 0; j <= bases[1].isize( ) - adaptn; j++ )
          {    if ( s2.Contains( adapt20, j ) )
               {    have_ad = True;
                    break;    }    }
          if (have_ad) return 0;
          cout << "past fourth check" << endl;    }
     Scram(0);    }
