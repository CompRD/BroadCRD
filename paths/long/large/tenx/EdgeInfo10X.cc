///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Show the barcodes of reads aligned to some or all of the given edges.

#include "MainTools.h"
#include "VecUtilities.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/tenx/TenxDirs.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(E, "comma-separated list of edge ids");
     CommandArgument_Int_Doc(N, "1 or 2 or 3 or 4 or 5 or 6");
     CommandArgument_Int_OrDefault_Doc(MODE, 1,
          "1 or 2; 1: output barcode; 2: output barcode.edge.edgepos");
     CommandArgument_Bool_OrDefault_Doc(AND, False,
          "if True, only show barcodes aligning to all of the edges");
     CommandArgument_Bool_OrDefault_Doc(POS_REL, False,
          "instead of reporting read positions, report position of read midpoint "
          "as fraction of edge length");
     EndCommandArguments;

     // Hardcoded constant: WARNING!!

     const int readlen = 88;

     // Hardcoded directories.

     String dir, odir, tdir;
     SetTenxDirs( N, dir, odir, tdir );

     // Parse edge ids.

     vec<int> es;
     ParseIntSet( "{" + E + "}", es );

     // Load data.

     HyperBasevectorX hb;
     BinaryReader::readFile( dir + "/a.hbx", &hb );
     vec<int> inv;
     BinaryReader::readFile( dir + "/a.inv", &inv );
     vec<vec<int>> ehits;
     vec<uint32_t> bcs;
     BinaryReader::readFile( odir + "/10X.ehits", &ehits );
     BinaryReader::readFile( tdir + "/10X.bc", &bcs );
     vec< pair<int,int> > places;
     if ( MODE == 2 ) BinaryReader::readFile( odir + "/10X.aligns", &places );

     // Find barcodes.

     vec<uint32_t> bx;
     vec<vec<uint32_t>> b( es.size( ) );
     for ( int j = 0; j < es.isize( ); j++ )
     {    int e = es[j];
          if ( e >= ehits.isize( ) )
          {    cout << "\nSomething is amiss.  Edge id too large." << endl;
               Scram(1);    }
          for ( int i = 0; i < ehits[e].isize( ); i++ )
               b[j].push_back( bcs[ ehits[e][i] ] );
          for ( int i = 0; i < ehits[ inv[e] ].isize( ); i++ )
               b[j].push_back( bcs[ ehits[ inv[e] ][i] ] );
          UniqueSort( b[j] );
          if ( !AND ) bx.append( b[j] );    }
     if (AND) Intersection( b, bx );
     else UniqueSort(bx);
     if ( MODE == 1 ) cout << printSeq(bx) << endl;
     if ( MODE == 2 )
     {    vec< triple<uint32_t,int,int> > x;
          for ( int j = 0; j < es.isize( ); j++ )
          {    int e = es[j];
               for ( int pass = 1; pass <= 2; pass++ )
               {    if ( pass == 2 ) e = inv[e];
                    for ( int i = 0; i < ehits[e].isize( ); i++ )
                    {    int rid = ehits[e][i];
                         uint32_t bc = bcs[rid];
                         if ( AND && !BinMember( bx, bc ) ) continue;
                         int pos = places[rid].second;
                         if ( pass == 2 ) 
                              pos = Max( 0, hb.Bases(e) - pos - readlen );
                         int ex = ( pass == 1 ? e : inv[e] );
                         x.push( bc, ex, pos );    }    }    }
          Sort(x);
          for ( int i = 0; i < x.isize( ); i++ )
          {    uint32_t bc = x[i].first;
               if ( i > 0 && bc != x[i-1].first ) cout << "\n";
               int e = x[i].second;
               cout << bc << ":" << e << ":";
               if ( !POS_REL ) cout << x[i].third;
               else
               {    double f = double( x[i].third + readlen/2 ) / hb.Bases(e);
                    if ( f < 0 ) f = 0;
                    if ( f > 1 ) f = 1;
                    cout << setprecision(3) << f;    }
               cout << "\n";    }
          cout << endl;    }
     Scram(0);    }
