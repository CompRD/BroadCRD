///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Align reads from one 10X barcode to a set of edges (normally comprising a single
// path).  Output is like EdgeInfo10X MODE=2 POS_REL=True.
//
// Mirrors Align10X.cc.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "kmers/KmerRecord.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/tenx/TenxDirs.h"
#include "paths/long/large/tenx/TenxTools.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_Int_Doc(N, "1 or 2 or 3 or 4 or 5 or 6");
     CommandArgument_Int_OrDefault_Doc(BC, -1, "barcode");
     CommandArgument_Int_OrDefault_Doc(RID, -1, "read id");
     CommandArgument_String_Doc(E, "comma-separated list of edges");
     EndCommandArguments;

     // Hardcoded readlength!!

     const int readlen = 88;

     // Check args.

     ForceAssert( ( BC >= 0 ) ^ ( RID >= 0 ) );

     // Hardcoded directories.

     String dir, odir, tdir;
     SetTenxDirs( N, dir, odir, tdir );

     // Parse edges.

     vec<int> es;
     ParseIntSet( E, es );

     // Load assembly.

     cout << Date( ) << ": loading assembly" << endl;
     HyperBasevector hb;
     BinaryReader::readFile( dir + "/a.hbv", &hb );
     vecbasevector tigs( dir + "/a.fastb" );
     vec<int> inv;
     BinaryReader::readFile( dir + "/a.inv", &inv );

     // Restrict to given edges.

     vec<int> esplus(es);
     for ( int i = 0; i < es.isize( ); i++ )
          esplus.push_back( inv[ es[i] ] );
     vec<int> dels;
     for ( int e = 0; e < (int) tigs.size( ); e++ )
     {    if ( !Member( esplus, e ) ) 
          {    tigs[e].resize(0);
               dels.push_back(e);    }    }
     hb.DeleteEdges(dels);

     // Read barcodes.

     vec<uint32_t> bcs;
     BinaryReader::readFile( tdir + "/10X.bc", &bcs );

     // Get barcode id.

     if ( RID >= 0 ) BC = bcs[RID];

     // Load 10X data.

     cout << Date( ) << ": loading 10X data" << endl;
     vec<int64_t> bci;
     BinaryReader::readFile( tdir + "/10X.bci", &bci );
     vecbasevector bases;
     bases.ReadRange( tdir + "/10X.fastb", bci[BC], bci[BC+1] );

     // Heuristics.

     const int max_freq = 5;
     const int max_frag = 1000;
     vec<int> KS     = {32,40,60,80,40,40};
     vec<int> starts = {0, 0, 0, 0, 40,10};
     const int Kbig = 200;

     // Create alignments.

     vec< pair<int,int> > places( bases.size( ), make_pair( -1, -1 ) );
     for ( int pass = 0; pass < KS.isize( ); pass++ )
     {    int K = KS[pass];
          if ( K == 32 ) 
          {    Match<32>( starts[pass], 
                    hb, tigs, inv, bases, Kbig, max_freq, max_frag, places );    }
          else if ( K == 40 ) 
          {    Match<40>( starts[pass], 
                    hb, tigs, inv, bases, Kbig, max_freq, max_frag, places );    }
          else if ( K == 60 ) 
          {    Match<60>( starts[pass], 
                    hb, tigs, inv, bases, Kbig, max_freq, max_frag, places );    }
          else if ( K == 80 ) 
          {    Match<80>( starts[pass], 
                    hb, tigs, inv, bases, Kbig, max_freq, max_frag, places );    }
          else
          {    cout << "Illegal K." << endl;
               Scram(1);    }    }

     // Report alignments.

     vec< quad<uint32_t,int,int,int> > x;
     for ( int i = 0; i < (int) bases.size( ); i++ )
     {    int rid = i + bci[BC];
          int e = places[i].first;
          if ( e < 0 ) continue;
          int pos = places[i].second;
          if ( Member( es, e ) ) x.push( BC, e, pos, rid );
          else x.push( BC, inv[e], hb.Bases(e) - pos - readlen, rid );    }
     Sort(x);
     for ( int i = 0; i < x.isize( ); i++ )
     {    uint32_t bc = x[i].first;
          if ( i > 0 && bc != x[i-1].first ) cout << "\n";
          int e = x[i].second;
          cout << bc << ":" << e << ":";
          double f = double( x[i].third + readlen/2 ) / hb.Bases(e);
          if ( f < 0 ) f = 0;
          if ( f > 1 ) f = 1;
          cout << setprecision(3) << f << " " << x[i].fourth << "\n";    }
     cout << endl;
     Scram(0);    }
