///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Find neighbors of a given line.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
// MakeDepend: library GMPXX

// require <cstddef> to use gmp in GCC-4.9
#include <cstddef>
#include <gmpxx.h>

#include "FastIfstream.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/Lines.h"
#include "paths/long/large/tenx/TenxDirs.h"
#include "paths/long/large/tenx/TenxTools.h"
#include <unordered_set>

typedef mpf_class big_float;

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_Int_Doc(S, "seed edge");
     CommandArgument_IntSet_Doc(E, "other edges");
     CommandArgument_Int_Doc(N, "dataset 1 or 2 or 3 or 4 or 5 or 6");
     EndCommandArguments;

     // Hardcoded directories.

     String dir, odir, tdir;
     SetTenxDirs( N, dir, odir, tdir );

     // Load counts and assembly.

     cout << "\n" << Date( ) << ": loading assembly" << endl;
     HyperBasevectorX hb;
     vec<int> inv, tol, lens;
     vec<vec<vec<vec<int>>>> lines;
     vec< vec< pair<int,int> > > aligns;
     #pragma omp parallel sections
     {    
          #pragma omp section
          {    BinaryReader::readFile( dir + "/a.hbx", &hb );    }
          #pragma omp section
          {    BinaryReader::readFile( dir + "/a.inv", &inv );    }
          #pragma omp section
          {    BinaryReader::readFile( dir + "/a.lines", &lines );    }
          #pragma omp section
          {    BinaryReader::readFile( dir + "/a.aligns", &aligns );    }    }

     // Load alignments, barcodes, and barcode indices

     cout << Date( ) << ": loading alignments and barcodes" << endl;
     vec< pair<int,int> > places;
     BinaryReader::readFile( odir + "/10X.aligns", &places );

     vec< uint32_t > bcs;
     BinaryReader::readFile( tdir + "/10X.bc", &bcs );

     vec< int64_t > bci;
     BinaryReader::readFile( tdir + "/10X.bci", &bci );


     cout << places.size() << " alignments, " << bcs.size()
             << " read barcodes, " << bci.size() << " unique barcodes" << endl;

     // figure out which barcodes impinge upon the seed edge
     unordered_set< uint32_t > targets;
     cout << Date( ) << ": finding barcodes on seed" << endl;
#pragma omp parallel for
     for ( size_t i = 0; i < places.size(); ++i ) {
         if ( places[i].first != -1  &&
                 ( places[i].first == S || places[i].first == inv[S] ) )
             targets.insert( bcs[i] );
     }
     cout << "there are " << targets.size() << " barcodes on the seed edge" << endl;


     vec<size_t> counts( E.size(), 0 );
     vec<vec<size_t>> hits_per_bc( E.size() );
     for ( auto& h : hits_per_bc ) h.resize(targets.size(),0);
     vec<uint32_t> targets1( targets.begin(), targets.end() );

     for ( size_t t = 0; t < targets1.size(); ++t ) {
         auto const target = targets1[t];
         // for each read with that barcode
         int64_t start = bci[target];
         ForceAssertGe(start, 0);
         int64_t end = ( target+1 < bci.size() ) ? bci[target+1] : bcs.size();
#pragma omp parallel for
         for ( int64_t i = start; i < end; ++i ) {
             // for each edge
             for ( size_t j = 0; j < E.size(); ++j )
                 if ( places[i].first != -1 && ( places[i].first == E[j] ||
                         places[i].first == inv[E[j]] )  ) {
#pragma omp atomic
                     counts[j]++;
#pragma omp atomic
                     hits_per_bc[j][t]++;
                 }
         }
     }

     for ( size_t j = 0; j < E.size(); ++j ) {
         cout << "edge " << E[j] << ": " << counts[j] << " hits: ";
         for ( size_t t = 0; t < targets.size(); ++t )
             cout << hits_per_bc[j][t] << " ";
         cout << endl;
     }


#if 0
     vec<int> es;
     ParseIntSet( E, es );

     vec<size_t> counts(es.size(), 0);

#pragma omp for
     for ( size_t i = 0; i < places.size(); ++i ) {
         for ( size_t j = 0; j < es.size(); ++j ) {
             if ( places[i].first == es[j] ||
                     places[i].first == inv[es[j]] )
#pragma omp atomic
                 counts[j]++;
         }
     }

     for ( size_t j = 0; j < es.size(); ++j ) {
         cout << "Edge " << es[j] << ", hits " << counts[j] << endl;
     }
#endif

     cout << "\n" << Date( ) << ": done\n" << endl;
     return 0;
}
