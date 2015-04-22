///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Show examples of magic E. coli reads for some random regions and their consensus.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "PrintAlignment.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "paths/long/magic/BasicScore.h"
#include "paths/long/magic/IterCon.h"
#include "paths/long/ultra/ThreadedBlocks.h"
#include "random/Random.h"

// QueryLookupTable K=12 MM=12 MC=0.005 SEQS=Jan30.2014.3/all.fastb L=scs.lookup 
//      SMITH_WAT=True FW_ONLY=True PARSEABLE=True > Jan30.2014.3/all.aligns

int main( )
{    RunTime( );

     SupportedHyperBasevector shb;
     BinaryReader::readFile( "/wga/scr4/macro/scs.shbv", &shb );

     vecbasevector magic;
     magic.ReadAll( "/wga/scr4/macro/Jan30.2014.3/all.fastb" );
     int N = magic.size( );

     vec<look_align> aligns;
     vec< vec<int> > aligns_index;
     LoadLookAligns( "/wga/scr4/macro/Jan30.2014.3/all.aligns", aligns,
          aligns_index, magic.size( ) );

     int len = 70;

     const int sample = 10;

     for ( int pass = 0; pass < sample; pass++ )
     {    int u;
          while(1)
          {    u = randomx( ) % shb.EdgeObjectCount( );
               if ( shb.EdgeLengthBases(u) >= 5000 ) break;    }
          int start = 1000 + ( randomx( ) % ( shb.EdgeLengthBases(u) - 2000 ) );
          int stop = start + len;
          cout << "\n==============================================================="
               << "=========================\n";
          cout << "\npass " << pass+1 << ", reference = " << len << "-base region on "
               << "DISCOVAR assembly edge " << u
               << ", starting at base " << start << endl;

          basevector b( shb.EdgeObject(u), start, len );
          basevector borig(b);

          vecbasevector segs;
          for ( int i = 0; i < aligns.isize( ); i++ )
          {    const look_align& la = aligns[i];
               int id1 = la.query_id;
               if ( aligns[i].target_id == u )
               {    if ( la.pos2( ) > start ) continue;
                    if ( la.Pos2( ) < stop ) continue;
                    int best_loc;
                    alignment a;
                    SmithWatFree( b, magic[id1], best_loc, a );
                    segs.push_back( 
                         basevector( magic[id1], a.pos2( ), a.Pos2( ) - a.pos2( ) ) );
                    align al(a);
                    cout << "\nalignment of read " << id1 << " to reference" << endl;
                    al.Flip( );
                    PrintVisualAlignment( False, cout, magic[id1], b, al );    }
               if ( aligns[i].target_id == shb.Inv(u) )
               {    align ar = la.a;
                    ar.ReverseThis( magic[id1].size( ), shb.EdgeLengthBases(u) );
                    if ( ar.pos2( ) > start ) continue;
                    if ( ar.Pos2( ) < stop ) continue;
                    int best_loc;
                    alignment a;
                    basevector mr = magic[id1];
                    mr.ReverseComplement( );
                    SmithWatFree( b, mr, best_loc, a );
                    segs.push_back( basevector( mr, a.pos2( ), a.Pos2( ) - a.pos2( ) ) );
                    cout << "\nalignment of read " << id1 << "' to reference" << endl;
                    align al(a);
                    al.Flip( );
                    PrintVisualAlignment( False, cout, mr, b, al );    }    }

          basevector con = IterCon(segs);
     
          int best_loc;
          alignment a;
          if ( con.size( ) <= b.size( ) )
               SmithWatFree( con, b, best_loc, a );
          else
          {    SmithWatFree( b, con, best_loc, a );
               a.Flip( );    }
          cout << "\nalignment of consensus to reference:\n";
          PrintVisualAlignment( False, cout, con, b, a );    }    }
