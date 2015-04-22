///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Align 10X reads.  Run after Bam10X.  For N=1, took ~30 minutes on crd26, would 
// run on half terabyte box.  Places 75.3%.
//
// Note that this keeps read ids as ints, could overflow.
//
// Generate:
// 1. 10X.aligns: vec<(edge,start on edge)> with one entry per read, set to -1
//                if no alignment was found
// 2. 10X.ehits:  for each edge, list of reads aligned to it.

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
     CommandArgument_String_OrDefault_Doc(HEAD, "10X", "head for output files");
     CommandArgument_Int_Doc(N, "1 or 2 or 3 or 4 or 5 or 6");
     EndCommandArguments;

     // Hardcoded directories.

     String dir, odir, tdir;
     SetTenxDirs( N, dir, odir, tdir );

     // Load assembly.

     cout << Date( ) << ": loading assembly" << endl;
     HyperBasevectorX hb;
     BinaryReader::readFile( dir + "/a.hbx", &hb );
     vecbasevector tigs( dir + "/a.fastb" );
     vec<int> inv;
     BinaryReader::readFile( dir + "/a.inv", &inv );

     // Load 10X data.

     cout << Date( ) << ": loading 10X data" << endl;
     vecbasevector bases( tdir + "/10X.fastb" );

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

     // Report stats and write alignments.

     int64_t count = 0;
     for ( int64_t i = 0; i < places.jsize( ); i++ )
          if ( places[i].first >= 0 ) count++;
     cout << PERCENT_RATIO( 3, count, places.jsize( ) ) << " of reads placed"
          << endl;
     cout << Date( ) << ": writing" << endl;
     BinaryWriter::writeFile( odir + "/" + HEAD + ".aligns", places );

     // Make edge hits.

     cout << Date( ) << ": generating edge hits" << endl;
     vec<vec<int>> ehits( tigs.size( ) );
     for ( int64_t i = 0; i < places.jsize( ); i++ )
          if ( places[i].first >= 0 ) ehits[ places[i].first ].push_back(i);
     BinaryWriter::writeFile( odir + "/" + HEAD + ".ehits", ehits );

     // Done.

     cout << Date( ) << ": done" << endl;
     cout << "peak memory usage = " << ToStringAddCommas( PeakMemUsageBytes( ) )
          << endl;
     Scram(0);    }
